# GitHub 与 Zenodo 上传指南

## 一、上传到 GitHub

### 1. 在 GitHub 建立空仓库

1. 登录 https://github.com/。
2. 点击右上角 `+`，选择 `New repository`。
3. Repository name 填写 `human-fascia-muscle-mechanosensitivity`。
4. 建议先选择 `Public`。如果投稿前暂时不公开，可先选择 `Private`，但最终 DOI 发布前应确认仓库访问策略。
5. 不要勾选自动添加 README、`.gitignore` 或 License，因为发布包中已经包含这些文件。
6. 点击 `Create repository`。

### 2. 从本地发布包推送

在 Windows Terminal、PowerShell 或 Git Bash 中进入发布目录：

```powershell
cd "PATH_TO_RELEASE_FOLDER"
git init
git config core.longpaths true
git add .
git commit -m "Release v1.0.0"
git branch -M main
git remote add origin https://github.com/YOUR_GITHUB_USERNAME/human-fascia-muscle-mechanosensitivity.git
git push -u origin main
```

将 `PATH_TO_RELEASE_FOLDER` 替换为本机发布目录，将 `YOUR_GITHUB_USERNAME` 替换为自己的 GitHub 用户名。首次推送时，GitHub 可能要求通过浏览器登录或使用 personal access token。不要把 token 写入脚本或仓库文件。

`git config core.longpaths true` 仅作用于当前仓库，用于避免 Windows 对较长结果路径报 `Filename too long`。不要省略这一行。

### 3. 建立正式 GitHub Release

1. 打开仓库主页，点击右侧 `Releases`。
2. 点击 `Draft a new release`。
3. 点击 `Choose a tag`，输入 `v1.0.0`，选择创建新标签。
4. Release title 填写 `Version 1.0.0`。
5. Release notes 可复制 `docs/RELEASE_NOTES.md` 的内容。
6. 点击 `Publish release`。

GitHub 会自动提供源代码 ZIP 和 TAR.GZ。也可以把项目生成的独立 ZIP 作为额外附件上传。

## 二、通过 GitHub–Zenodo 集成归档

### 1. 连接账户

1. 登录 https://zenodo.org/，建议使用 ORCID 或 GitHub 登录。
2. 打开账户设置中的 GitHub integration 页面。
3. 授权 Zenodo 访问 GitHub。
4. 在仓库列表中找到 `human-fascia-muscle-mechanosensitivity` 并开启归档开关。

应在发布 GitHub Release 之前完成这一步。若已经先发布 Release，可在连接后重新创建一个补丁版本，例如 `v1.0.1`，或按 Zenodo 页面提示重新触发归档。

### 2. 触发 Zenodo 归档

在 GitHub 发布 `v1.0.0` Release 后，Zenodo 会自动建立记录草稿或存档。进入 Zenodo 的 Uploads 页面检查：

- Title：`Human fascia–muscle mechanosensitivity transcriptomic analysis`
- Upload type：`Software`
- Version：`1.0.0`
- Creators：Cheng Chen；Hai Huang；Jian Sun
- Jian Sun ORCID：`0000-0001-7189-8493`
- Description、keywords 和 licence：参考 `docs/ZENODO_METADATA.md`
- Related identifier：填写 GitHub 仓库 URL

确认文件、作者顺序和描述后发布 Zenodo 记录。发布后 Zenodo 会生成：

- version DOI：对应 `v1.0.0`；
- concept DOI：对应整个项目的全部未来版本。

论文中通常引用本次分析所对应的 version DOI。

## 三、如果不使用 GitHub 集成

也可以在 Zenodo 点击 `New upload`，直接上传本目录旁边生成的 ZIP 包。然后按 `docs/ZENODO_METADATA.md` 填写元数据并发布。采用这种方式时，要在 Related identifiers 中手动加入 GitHub 仓库地址。

## 四、取得 DOI 后更新论文

将 Data Availability 中的未来时态替换为：

```text
The analysis code, frozen registries, source inventories, analysis decision records and machine-readable result tables are available in a version-controlled GitHub repository at [GITHUB URL]. Version 1.0.0 has been archived in Zenodo at https://doi.org/[ZENODO DOI]. The analysis code is released under the MIT License.
```

把 `[GITHUB URL]` 和 `[ZENODO DOI]` 替换为真实地址。不要在 Zenodo 正式分配 DOI 之前自行填写或猜测 DOI。

## 五、最终核验

1. 在未登录 GitHub 的浏览器窗口中打开仓库，确认公开可见。
2. 点击 Zenodo DOI，确认能打开正确记录。
3. 下载 Zenodo ZIP，核对 README、脚本、结果表和 6 张图是否存在。
4. 用 `metadata/SHA256SUMS.txt` 抽查下载文件。
5. 确认仓库中没有原始 RDS/H5AD、缓存目录、账号 token、个人绝对路径或 Word 投稿稿件。
