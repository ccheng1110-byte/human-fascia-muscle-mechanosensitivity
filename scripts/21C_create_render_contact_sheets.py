"""Create four-page contact sheets for visual QA of the rendered manuscript."""

from pathlib import Path

from PIL import Image, ImageDraw


render_dir = Path(
    r"./results/11_manuscript_preparation/"
    r"21C_submission_strengthened_word_manuscript/rendered_QA"
)

for start in (1, 5, 9, 13):
    canvas = Image.new("RGB", (1440, 2040), "white")
    draw = ImageDraw.Draw(canvas)
    for index, page_number in enumerate(range(start, start + 4)):
        page = Image.open(render_dir / f"page-{page_number}.png").convert("RGB")
        page.thumbnail((700, 950))
        x = (index % 2) * 720 + (720 - page.width) // 2
        y = (index // 2) * 1020 + 45
        canvas.paste(page, (x, y))
        draw.text((x, 15 + (index // 2) * 1020), f"Page {page_number}", fill="black")
    canvas.save(render_dir / f"contact-{start:02d}-{start + 3:02d}.png")
