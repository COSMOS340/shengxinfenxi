#!/usr/bin/env python3
import argparse
import ctypes
import datetime
import json
import os
import platform
import shutil
import socket
import subprocess


class MemoryStatusEx(ctypes.Structure):
    _fields_ = [
        ("dwLength", ctypes.c_ulong),
        ("dwMemoryLoad", ctypes.c_ulong),
        ("ullTotalPhys", ctypes.c_ulonglong),
        ("ullAvailPhys", ctypes.c_ulonglong),
        ("ullTotalPageFile", ctypes.c_ulonglong),
        ("ullAvailPageFile", ctypes.c_ulonglong),
        ("ullTotalVirtual", ctypes.c_ulonglong),
        ("ullAvailVirtual", ctypes.c_ulonglong),
        ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
    ]


def windows_memory():
    status = MemoryStatusEx()
    status.dwLength = ctypes.sizeof(MemoryStatusEx)
    ok = ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status))
    if not ok:
        raise ctypes.WinError()
    return {
        "total_physical_bytes": status.ullTotalPhys,
        "available_physical_bytes": status.ullAvailPhys,
        "memory_load_percent": status.dwMemoryLoad,
    }


def gpu_snapshot():
    executable = shutil.which("nvidia-smi")
    if executable is None:
        return {"available": False, "reason": "nvidia-smi not found"}
    command = [
        executable,
        "--query-gpu=name,memory.total,memory.free,driver_version",
        "--format=csv,noheader,nounits",
    ]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        return {
            "available": False,
            "reason": completed.stderr.strip() or "nvidia-smi failed",
        }
    rows = []
    for line in completed.stdout.splitlines():
        values = [value.strip() for value in line.split(",")]
        if len(values) == 4:
            rows.append({
                "name": values[0],
                "memory_total_mib": values[1],
                "memory_free_mib": values[2],
                "driver_version": values[3],
            })
    return {"available": True, "devices": rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--workdir", required=True)
    args = parser.parse_args()

    disk = shutil.disk_usage(args.workdir)
    payload = {
        "captured_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "logical_cpu_count": os.cpu_count(),
        "memory": windows_memory(),
        "project_disk": {
            "path": os.path.abspath(args.workdir),
            "total_bytes": disk.total,
            "used_bytes": disk.used,
            "free_bytes": disk.free,
        },
        "gpu": gpu_snapshot(),
    }
    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
