"""Export docs to GitHub Wiki Markdown. Does not access GitHub or push changes."""
import argparse
from pathlib import Path
import re
import shutil
from urllib.parse import urlparse


def export(repo_url, destination):
    root = Path(__file__).resolve().parents[1]
    destination = Path(destination).resolve()
    parsed = urlparse(repo_url.rstrip("/"))
    parts = parsed.path.strip("/").split("/")
    if parsed.scheme != "https" or parsed.netloc != "github.com" or len(parts) != 2:
        raise ValueError("Use https://github.com/USERNAME/REPOSITORY")
    repo_url = repo_url.rstrip("/").removesuffix(".git")
    if destination == root or destination == root / "docs":
        raise ValueError("Destination must not be the repository root or docs directory")
    destination.mkdir(parents=True, exist_ok=True)
    pages = {p.stem for p in (root / "docs").glob("*.md") if p.stem != "Validation"}
    def convert(match):
        label, link = match.group(1), match.group(2)
        if link.startswith(("https://", "http://", "#")):
            return match.group(0)
        path, _, anchor = link.partition("#")
        fragment = "#" + anchor if anchor else ""
        if path.startswith("assets/"):
            target = repo_url + "/raw/refs/heads/main/docs/" + path
        elif path == "../README.md":
            target = "Home"
        elif path.endswith(".md") and Path(path).stem in pages and not path.startswith("../"):
            target = Path(path).stem
        else:
            target = repo_url + "/blob/main/" + (path[3:] if path.startswith("../") else "docs/" + path)
        return f"[{label}]({target}{fragment})"
    for page in (root / "docs").glob("*.md"):
        if page.stem == "Validation":
            continue
        content = re.sub(r"\[([^\]]*)\]\(([^)]+)\)", convert, page.read_text())
        (destination / page.name).write_text(content)
    home = """# K-tools Wiki

K-tools extracts RNA k-mer signatures and associates them with RBP profiles.

- [Documentation](Documentation)
- [Installation](Installation)
- [Tutorial: simulated RNA sequences](Tutorial)
- [miRNA seed annotation](MiRNA-seeds)
- [CDR1as example](CDR1as-miRNA-example)
- [Kmap](Kmap)
- [CDS domain enrichment](Domain-enrichment)
- [Methodology](Methodology)
- [Detailed module parameters](Module-parameters)
- [API reference](API-reference)
- [Benchmarks: in silico and eCLIP](Benchmarks)
- [Data and reproducibility](Data-and-reproducibility)
- [FAQs](FAQs)
- [Known limitations](Known-limitations)
- [Citation](Citation)

[Source code and README](REPO_URL)
"""
    (destination / "Home.md").write_text(home.replace("REPO_URL", repo_url))
    (destination / "_Sidebar.md").write_text("[K-tools](Home)\n\n" + "\n".join(
        f"- [{name.replace('-', ' ')}]({name})" for name in [
            "Documentation", "Installation", "Tutorial", "Methodology", "API-reference",
            "Benchmarks", "Kmap", "Domain-enrichment", "MiRNA-seeds", "CDR1as-miRNA-example", "Module-parameters", "FAQs", "Data-and-reproducibility", "Known-limitations", "Citation"]) + "\n")
    (destination / "_Footer.md").write_text(f"[Repository]({repo_url}) · [Tutorial](Tutorial) · [FAQs](FAQs)\n")
    print(f"Wiki pages exported to {destination}; no network action performed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-url", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    export(args.repo_url, args.out)
