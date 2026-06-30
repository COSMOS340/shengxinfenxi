#!/usr/bin/env python3
"""Create a bioinformatics reproduction project scaffold from a route matrix."""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_VERSION = "2026-06-30"

ROUTE_COLUMNS = [
    "route_id",
    "route_name",
    "required_data_dirs",
    "required_metadata",
    "required_logs",
    "required_method_modules",
    "figure_playbook",
]

DIR_ROOT_MAP = {
    "metadata": "00_metadata",
    "plan": "01_plan",
    "scripts": "02_scripts",
    "data_raw": "03_data_raw",
    "data_processed": "04_data_processed",
    "results": "05_results",
    "figures": "06_figures",
    "tables": "07_tables",
    "logs": "99_logs",
}

CORE_DIRS = [
    "00_metadata",
    "01_plan",
    "02_scripts",
    "02_scripts/00_setup",
    "02_scripts/01_data_download",
    "02_scripts/02_preprocessing",
    "02_scripts/03_discovery",
    "02_scripts/04_modeling",
    "02_scripts/05_validation",
    "02_scripts/06_mechanism",
    "02_scripts/07_figures",
    "03_data_raw",
    "04_data_processed",
    "05_results",
    "06_figures",
    "06_figures/main",
    "06_figures/supplementary",
    "06_figures/review",
    "07_tables",
    "07_tables/main",
    "07_tables/supplementary",
    "99_logs",
]

TEMPLATE_DESTINATIONS = {
    "analysis_plan.md": "01_plan/analysis_plan.md",
    "limitations.md": "01_plan/limitations.md",
    "data_manifest.tsv": "00_metadata/data_manifest.tsv",
    "sample_metadata.tsv": "00_metadata/sample_metadata.tsv",
    "methods_log.tsv": "99_logs/methods_log.tsv",
    "figure_manifest.tsv": "99_logs/figure_manifest.tsv",
    "validation_log.tsv": "99_logs/validation_log.tsv",
    "handoff_checklist.md": "99_logs/handoff_checklist.md",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a route-specific bioinformatics reproduction scaffold."
    )
    parser.add_argument("--route-id", help="Exact route_id from the route matrix, for example R02.")
    parser.add_argument("--project-name", help="Name of the project directory to create.")
    parser.add_argument(
        "--output-root",
        default=".",
        help="Directory where the project directory will be created. Defaults to the current directory.",
    )
    parser.add_argument("--route-matrix", help="Path to route_scaffold_matrix.tsv.")
    parser.add_argument("--template-dir", help="Path to reproduction scaffold templates.")
    parser.add_argument("--list-routes", action="store_true", help="Print available routes and exit.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned actions without writing files.")
    return parser.parse_args()


def split_cells(value: str) -> list[str]:
    return [part.strip() for part in value.split(";") if part.strip()]


def script_dir() -> Path:
    return Path(__file__).resolve().parent


def default_route_matrix_paths() -> list[Path]:
    base = script_dir()
    return [
        base.parent.parent / "99_logs" / "reproduction_project_scaffold" / "route_scaffold_matrix.tsv",
        base.parent / "assets" / "reproduction-project-scaffold" / "route_scaffold_matrix.tsv",
        base.parent.parent / "assets" / "reproduction-project-scaffold" / "route_scaffold_matrix.tsv",
    ]


def default_template_paths() -> list[Path]:
    base = script_dir()
    return [
        base.parent / "templates" / "reproduction_project_scaffold",
        base.parent / "assets" / "reproduction-project-scaffold",
        base.parent.parent / "assets" / "reproduction-project-scaffold",
    ]


def resolve_path(user_value: str | None, search_paths: list[Path], label: str) -> Path:
    if user_value:
        path = Path(user_value).expanduser().resolve()
        if not path.exists():
            raise SystemExit(f"{label} does not exist: {path}")
        return path

    for path in search_paths:
        if path.exists():
            return path.resolve()

    checked = "\n".join(str(path) for path in search_paths)
    raise SystemExit(f"{label} was not found. Checked:\n{checked}")


def load_routes(route_matrix: Path) -> dict[str, dict[str, str]]:
    with route_matrix.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"Route matrix is empty: {route_matrix}")
        missing = [column for column in ROUTE_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise SystemExit(f"Route matrix lacks columns: {', '.join(missing)}")
        routes = {row["route_id"]: row for row in reader if row.get("route_id")}

    if not routes:
        raise SystemExit(f"Route matrix has no route rows: {route_matrix}")
    return routes


def map_logical_dir(logical_dir: str) -> str:
    parts = [part for part in logical_dir.strip("/").split("/") if part]
    if not parts:
        raise SystemExit("Route matrix contains an empty directory value.")
    root = parts[0]
    if root not in DIR_ROOT_MAP:
        raise SystemExit(f"Route directory uses an unsupported root: {logical_dir}")
    return "/".join([DIR_ROOT_MAP[root], *parts[1:]])


def validate_project_name(project_name: str) -> None:
    if project_name in {"", ".", ".."}:
        raise SystemExit("Project name must be a directory name.")
    path = Path(project_name)
    if path.name != project_name or path.is_absolute():
        raise SystemExit("Project name must not contain path separators.")


def ensure_empty_target(target: Path) -> None:
    if target.exists() and any(target.iterdir()):
        raise SystemExit(f"Target directory exists and is not empty: {target}")


def require_templates(template_dir: Path) -> None:
    missing = [name for name in TEMPLATE_DESTINATIONS if not (template_dir / name).is_file()]
    if missing:
        raise SystemExit(f"Template directory lacks files: {', '.join(missing)}")


def write_text(path: Path, content: str, dry_run: bool) -> None:
    if dry_run:
        print(f"write {path}")
        return
    path.write_text(content, encoding="utf-8")


def make_dir(path: Path, dry_run: bool) -> None:
    if dry_run:
        print(f"mkdir {path}")
        return
    path.mkdir(parents=True, exist_ok=True)


def copy_file(source: Path, dest: Path, dry_run: bool) -> None:
    if dry_run:
        print(f"copy {source} -> {dest}")
        return
    shutil.copy2(source, dest)


def make_empty_file(path: Path, dry_run: bool) -> None:
    if dry_run:
        print(f"touch {path}")
        return
    path.touch(exist_ok=True)


def route_summary(project_name: str, route: dict[str, str], created_dirs: list[str]) -> str:
    method_modules = split_cells(route["required_method_modules"])
    metadata_files = split_cells(route["required_metadata"])
    log_files = split_cells(route["required_logs"])
    lines = [
        "# Route Summary",
        "",
        f"- Project name: {project_name}",
        f"- Route ID: {route['route_id']}",
        f"- Route name: {route['route_name']}",
        f"- Figure playbook: {route['figure_playbook']}",
        "",
        "## Required Method Modules",
        "",
    ]
    lines.extend(f"- {name}" for name in method_modules)
    lines.extend(["", "## Required Metadata Files", ""])
    lines.extend(f"- {name}" for name in metadata_files)
    lines.extend(["", "## Required Log Files", ""])
    lines.extend(f"- {name}" for name in log_files)
    lines.extend(["", "## Created Route Directories", ""])
    lines.extend(f"- {name}" for name in created_dirs)
    lines.append("")
    return "\n".join(lines)


def generation_manifest(
    project_name: str,
    route: dict[str, str],
    route_matrix: Path,
    template_dir: Path,
    created_dirs: list[str],
) -> str:
    rows = [
        ("script_version", SCRIPT_VERSION),
        ("generated_at", datetime.now(timezone.utc).isoformat()),
        ("project_name", project_name),
        ("route_id", route["route_id"]),
        ("route_name", route["route_name"]),
        ("route_matrix_path", str(route_matrix)),
        ("template_dir", str(template_dir)),
        ("required_data_dirs", route["required_data_dirs"]),
        ("created_data_dirs", ";".join(created_dirs)),
        ("required_metadata", route["required_metadata"]),
        ("required_logs", route["required_logs"]),
        ("required_method_modules", route["required_method_modules"]),
        ("figure_playbook", route["figure_playbook"]),
    ]
    return "field\tvalue\n" + "\n".join(f"{key}\t{value}" for key, value in rows) + "\n"


def create_project(args: argparse.Namespace, route: dict[str, str], route_matrix: Path, template_dir: Path) -> Path:
    project_name = args.project_name
    if project_name is None:
        raise SystemExit("--project-name is required unless --list-routes is used.")
    validate_project_name(project_name)

    output_root = Path(args.output_root).expanduser().resolve()
    target = output_root / project_name
    if not args.dry_run:
        output_root.mkdir(parents=True, exist_ok=True)
    ensure_empty_target(target)

    require_templates(template_dir)

    route_dirs = [map_logical_dir(value) for value in split_cells(route["required_data_dirs"])]
    all_dirs = list(dict.fromkeys([*CORE_DIRS, *route_dirs]))

    make_dir(target, args.dry_run)
    for rel_path in all_dirs:
        make_dir(target / rel_path, args.dry_run)

    for template_name, rel_dest in TEMPLATE_DESTINATIONS.items():
        copy_file(template_dir / template_name, target / rel_dest, args.dry_run)

    copied_names = {
        "data_manifest.tsv",
        "sample_metadata.tsv",
        "methods_log.tsv",
        "figure_manifest.tsv",
        "validation_log.tsv",
    }
    for name in split_cells(route["required_metadata"]):
        if name not in copied_names:
            make_empty_file(target / "00_metadata" / name, args.dry_run)
    for name in split_cells(route["required_logs"]):
        if name not in copied_names:
            make_empty_file(target / "99_logs" / name, args.dry_run)

    write_text(target / "01_plan" / "route_summary.md", route_summary(project_name, route, route_dirs), args.dry_run)
    write_text(
        target / "99_logs" / "scaffold_generation_manifest.tsv",
        generation_manifest(project_name, route, route_matrix, template_dir, route_dirs),
        args.dry_run,
    )

    return target


def main() -> int:
    args = parse_args()
    route_matrix = resolve_path(args.route_matrix, default_route_matrix_paths(), "Route matrix")
    template_dir = resolve_path(args.template_dir, default_template_paths(), "Template directory")
    routes = load_routes(route_matrix)

    if args.list_routes:
        print("route_id\troute_name")
        for route_id in sorted(routes):
            print(f"{route_id}\t{routes[route_id]['route_name']}")
        return 0

    if args.route_id is None:
        raise SystemExit("--route-id is required unless --list-routes is used.")
    if args.route_id not in routes:
        route_ids = ", ".join(sorted(routes))
        raise SystemExit(f"Route ID not found: {args.route_id}. Available route IDs: {route_ids}")

    target = create_project(args, routes[args.route_id], route_matrix, template_dir)
    if args.dry_run:
        print(f"dry_run_target\t{target}")
    else:
        print(f"created_project\t{target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
