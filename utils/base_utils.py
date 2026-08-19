import re

def get_current_version(changelog_path="CHANGELOG.md"):
    with open(changelog_path, "r", encoding="utf-8") as f:
        content = f.read()

    match = re.search(r"## \[(\d+\.\d+\.\d+)\]", content)
    if match:
        return match.group(1)

    return "unknown"