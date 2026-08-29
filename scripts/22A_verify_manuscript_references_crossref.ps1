$ErrorActionPreference = "Stop"

$project = "."
$manuscript = Join-Path $project "results\11_manuscript_preparation\22B_style_polished_manuscript_v7.md"
$outputDir = Join-Path $project "results\11_manuscript_preparation\22A_final_reference_verification"
$outputCsv = Join-Path $outputDir "22A_reference_field_verification_v2.csv"
$outputMd = Join-Path $outputDir "22A_reference_verification_summary_v2.md"

New-Item -ItemType Directory -Force -Path $outputDir | Out-Null

function Normalize-Text([string]$value) {
    if ([string]::IsNullOrWhiteSpace($value)) { return "" }
    $formD = $value.Normalize([Text.NormalizationForm]::FormD)
    $withoutMarks = [Regex]::Replace($formD, "\p{Mn}", "")
    return ([Regex]::Replace($withoutMarks.ToLowerInvariant(), "[^a-z0-9]+", "")).Trim()
}

function Get-FirstFamilies([string]$authors) {
    $authorText = ($authors -replace ",\s*et al\.?$", "").Trim()
    $parts = $authorText -split ","
    $families = foreach ($part in $parts) {
        $clean = ($part -replace "\s+", " ").Trim()
        if ($clean) {
            $tokens = @($clean -split " ")
            if ($tokens.Count -gt 1 -and $tokens[-1] -match "^[A-Z]{1,4}\.?$") {
                $family = ($tokens[0..($tokens.Count - 2)] -join " ")
            } else {
                $family = $tokens[0]
            }
            Normalize-Text $family
        }
    }
    return @($families | Where-Object { $_ })
}

function Get-CrossrefYear($message) {
    foreach ($field in @("published-print", "published-online", "published", "issued")) {
        if ($message.PSObject.Properties.Name -contains $field) {
            $parts = $message.$field.'date-parts'
            if ($parts -and $parts[0] -and $parts[0][0]) { return [int]$parts[0][0] }
        }
    }
    return $null
}

function Get-CrossrefValue($message, [string]$field) {
    if ($message.PSObject.Properties.Name -contains $field) {
        $value = $message.$field
        if ($value -is [array]) { return [string]$value[0] }
        return [string]$value
    }
    return ""
}

function Test-TitleMatch([string]$current, [string]$source) {
    $a = Normalize-Text $current
    $b = Normalize-Text $source
    if ($a -eq $b) { return $true }
    if (-not $a -or -not $b) { return $false }
    return ($a.Contains($b) -or $b.Contains($a))
}

function Test-JournalMatch([string]$current, [string]$source) {
    $a = Normalize-Text $current
    $b = Normalize-Text $source
    return ($a -eq $b -or $a.Contains($b) -or $b.Contains($a))
}

function Test-PageMatch([string]$current, [string]$sourcePage, [string]$sourceArticle) {
    $a = Normalize-Text $current
    $b = Normalize-Text $sourcePage
    $c = Normalize-Text $sourceArticle
    if (-not $a) { return $false }
    return ($a -eq $b -or $a -eq $c -or ($b -and $b.Contains($a)) -or ($c -and $c.Contains($a)))
}

$lines = Get-Content -LiteralPath $manuscript -Encoding UTF8
$referenceLines = $lines | Where-Object { $_ -match '^\d+\. ' }
$rows = @()

foreach ($line in $referenceLines) {
    $match = [Regex]::Match(
        $line,
        '^(?<id>\d+)\.\s+(?<authors>.*?)\.\s+(?<title>.*?)\.\s+\*(?<journal>.*?)\*\.\s+(?<year>\d{4});(?<volume>[^:]+):(?<pages>.*?)\.\s+doi:(?<doi>\S+)$'
    )
    if (-not $match.Success) {
        $rows += [pscustomobject]@{ id = "?"; doi = ""; http_status = "parse_error"; overall_status = "Needs manual check"; notes = "Reference line could not be parsed: $line" }
        continue
    }

    $id = [int]$match.Groups["id"].Value
    $doi = $match.Groups["doi"].Value.Trim().TrimEnd('.')
    $currentAuthors = $match.Groups["authors"].Value.Trim()
    $currentTitle = $match.Groups["title"].Value.Trim()
    $currentJournal = $match.Groups["journal"].Value.Trim()
    $currentYear = [int]$match.Groups["year"].Value
    $currentVolume = $match.Groups["volume"].Value.Trim()
    $currentPages = $match.Groups["pages"].Value.Trim()

    try {
        $response = Invoke-WebRequest -Uri ("https://api.crossref.org/works/" + $doi) -Headers @{ "User-Agent" = "manuscript-reference-audit/1.0" } -TimeoutSec 45 -UseBasicParsing
        $sourceMessage = ($response.Content | ConvertFrom-Json).message
    } catch {
        $rows += [pscustomobject]@{
            id = $id; doi = $doi; http_status = "HTTP_ERROR"; doi_resolves = $false
            current_first_author = (Get-FirstFamilies $currentAuthors)[0]; crossref_first_author = ""
            current_authors_first_three = ((Get-FirstFamilies $currentAuthors) -join ";"); crossref_authors_first_three = ""
            current_title = $currentTitle; crossref_title = ""; current_journal = $currentJournal; crossref_journal = ""
            current_year = $currentYear; crossref_year = ""; current_volume = $currentVolume; crossref_volume = ""
            current_issue = "not_reported"; crossref_issue = ""; current_pages_or_article = $currentPages; crossref_pages_or_article = ""
            title_match = $false; authors_match = $false; journal_match = $false; year_match = $false; volume_match = $false; pages_match = $false
            overall_status = "Needs fix or manual check"; notes = $_.Exception.Message
        }
        continue
    }

    $sourceAllFamilies = @($sourceMessage.author | ForEach-Object { Normalize-Text $_.family })
    $sourceFamilies = @($sourceAllFamilies | Select-Object -First 3)
    $currentFamilies = @(Get-FirstFamilies $currentAuthors)
    $sourceTitle = Get-CrossrefValue $sourceMessage "title"
    $sourceJournal = Get-CrossrefValue $sourceMessage "container-title"
    $sourceYear = Get-CrossrefYear $sourceMessage
    $sourceVolume = Get-CrossrefValue $sourceMessage "volume"
    $sourceIssue = Get-CrossrefValue $sourceMessage "issue"
    $sourcePage = Get-CrossrefValue $sourceMessage "page"
    $sourceArticle = Get-CrossrefValue $sourceMessage "article-number"
    $sourcePagesOrArticle = if ($sourcePage) { $sourcePage } else { $sourceArticle }

    $titleMatch = Test-TitleMatch $currentTitle $sourceTitle
    if ($currentAuthors -match '(?i)\bet al\.?') {
        $authorsToCompare = [array]($sourceAllFamilies | Select-Object -First $currentFamilies.Count)
    } else {
        $authorsToCompare = [array]$sourceAllFamilies
    }
    $authorsToCompare = @($authorsToCompare)
    $authorsMatch = ($currentFamilies.Count -gt 0 -and $authorsToCompare.Count -ge $currentFamilies.Count)
    if ($authorsMatch) {
        for ($i = 0; $i -lt $currentFamilies.Count; $i++) {
            if ($currentFamilies[$i] -ne $authorsToCompare[$i]) { $authorsMatch = $false; break }
        }
    }
    $journalMatch = Test-JournalMatch $currentJournal $sourceJournal
    $yearMatch = ($currentYear -eq $sourceYear)
    $volumeMatch = ($currentVolume -eq $sourceVolume)
    $pagesMatch = Test-PageMatch $currentPages $sourcePage $sourceArticle
    $allMatch = ($titleMatch -and $authorsMatch -and $journalMatch -and $yearMatch -and $volumeMatch -and $pagesMatch)
    $notes = @()
    if (-not $titleMatch) { $notes += "title" }
    if (-not $authorsMatch) { $notes += "author order/name" }
    if (-not $journalMatch) { $notes += "journal" }
    if (-not $yearMatch) { $notes += "year" }
    if (-not $volumeMatch) { $notes += "volume" }
    if (-not $pagesMatch) { $notes += "pages/article number" }
    $rows += [pscustomobject]@{
        id = $id; doi = $doi; http_status = 200; doi_resolves = $true
        current_first_author = $currentFamilies[0]; crossref_first_author = $sourceFamilies[0]
        current_authors_first_three = ($currentFamilies -join ";"); crossref_authors_first_three = ($sourceFamilies -join ";")
        current_title = $currentTitle; crossref_title = $sourceTitle
        current_journal = $currentJournal; crossref_journal = $sourceJournal
        current_year = $currentYear; crossref_year = $sourceYear
        current_volume = $currentVolume; crossref_volume = $sourceVolume
        current_issue = "not_reported"; crossref_issue = $sourceIssue
        current_pages_or_article = $currentPages; crossref_pages_or_article = $sourcePagesOrArticle
        title_match = $titleMatch; authors_match = $authorsMatch; journal_match = $journalMatch; year_match = $yearMatch; volume_match = $volumeMatch; pages_match = $pagesMatch
        overall_status = if ($allMatch) { "Verified" } else { "Check suggested" }
        notes = ($notes -join "; ")
    }
    Start-Sleep -Milliseconds 150
}

$rows | Export-Csv -LiteralPath $outputCsv -NoTypeInformation -Encoding UTF8
$verified = @($rows | Where-Object { $_.overall_status -eq "Verified" }).Count
$check = @($rows | Where-Object { $_.overall_status -eq "Check suggested" }).Count
$errors = @($rows | Where-Object { $_.overall_status -eq "Needs fix or manual check" }).Count
$doiResolved = @($rows | Where-Object { $_.doi_resolves -eq $true }).Count
$report = @(
    "# Final reference verification",
    "",
    "Verification route: Crossref DOI metadata, queried on $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss K').",
    "This verifies DOI resolvability and bibliographic metadata. It does not imply that every article is open access or that every full text is freely downloadable.",
    "",
    "- Total references checked: $($rows.Count)",
    "- DOI metadata records resolved: $doiResolved/$($rows.Count)",
    "- Verified across title, listed author prefix/order, journal, year, volume and pages/article number: $verified",
    "- Check suggested: $check",
    "- Needs fix or manual check: $errors",
    "",
    "## Interpretation",
    "",
    "A `Verified` row means the manuscript citation agrees with Crossref on the checked fields after normalization of punctuation, diacritics and journal-name variants. `Check suggested` rows require manual inspection, especially where Crossref has incomplete page or article-number metadata."
)
if ($check -gt 0 -or $errors -gt 0) {
    $report += ""
    $report += "## Rows requiring attention"
    $report += ""
    $report += "| ID | DOI | Status | Fields |"
    $report += "|---:|---|---|---|"
    foreach ($row in ($rows | Where-Object { $_.overall_status -ne "Verified" })) {
        $report += "| $($row.id) | $($row.doi) | $($row.overall_status) | $($row.notes) |"
    }
}
Set-Content -LiteralPath $outputMd -Value $report -Encoding UTF8
Write-Output "Final reference verification completed."
Write-Output "CSV: $outputCsv"
Write-Output "Summary: $outputMd"
Write-Output "Resolved: $doiResolved/$($rows.Count); verified: $verified; check suggested: $check; errors: $errors"
