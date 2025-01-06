# PFEIFER - MAPD - Modified by FrogAi for FrogPilot
#!/usr/bin/env python3
import json
import os
import stat
import subprocess
import time
import urllib.request

from pathlib import Path

from openpilot.selfdrive.frogpilot.frogpilot_utilities import is_url_pingable
from openpilot.selfdrive.frogpilot.frogpilot_variables import MAPD_PATH, MAPS_PATH

VERSION = "v1"

GITHUB_VERSION_URL = f"https://github.com/FrogAi/FrogPilot-Resources/raw/Versions/mapd_version_{VERSION}.json"
GITLAB_VERSION_URL = f"https://gitlab.com/FrogAi/FrogPilot-Resources/-/raw/Versions/mapd_version_{VERSION}.json"

VERSION_PATH = Path("/data/media/0/osm/mapd_version")

def is_mapd_running():
  try:
    result = subprocess.run(["pgrep", "-f", str(MAPD_PATH)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.returncode == 0
  except Exception as e:
    print(f"Error checking if mapd is running: {e}")
    return False

def download(current_version):
  urls = [
    f"https://github.com/pfeiferj/openpilot-mapd/releases/download/{current_version}/mapd",
    f"https://gitlab.com/FrogAi/FrogPilot-Resources/-/raw/Mapd/{current_version}"
  ]

  MAPD_PATH.parent.mkdir(parents=True, exist_ok=True)
  tmp_path = MAPD_PATH.with_suffix(".tmp")

  for url in urls:
    try:
      with urllib.request.urlopen(url, timeout=5) as f:
        data = f.read()
        with tmp_path.open('wb') as output:
          output.write(data)

      tmp_path.chmod(tmp_path.stat().st_mode | stat.S_IEXEC)

      os.replace(tmp_path, MAPD_PATH)
      VERSION_PATH.write_text(current_version)
      print(f"Updated mapd to version {current_version} from {url}")
      return
    except Exception as error:
      print(f"Failed to download mapd from {url}: {error}")

def get_installed_version():
  if not VERSION_PATH.exists():
    return "v0"
  return VERSION_PATH.read_text().strip()

def get_latest_version():
  for url in [GITHUB_VERSION_URL, GITLAB_VERSION_URL]:
    try:
      with urllib.request.urlopen(url, timeout=5) as response:
        return json.loads(response.read().decode('utf-8'))['version']
    except Exception as error:
      print(f"Failed to get version from {url}: {error}")
  return "v0"

def update_mapd():
  installed_version = get_installed_version()
  latest_version = get_latest_version()

  if latest_version == "v0":
    return

  if installed_version != latest_version:
    print("New mapd version available, stopping the mapd process for update")
    subprocess.run(["pkill", "-f", str(MAPD_PATH)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    download(latest_version)
  else:
    print("Mapd is up to date")

def ensure_mapd_is_running():
  if is_mapd_running():
    print("mapd is already running.")
    return

  if MAPD_PATH.exists():
    if not os.access(MAPD_PATH, os.X_OK):
      try:
        MAPD_PATH.chmod(MAPD_PATH.stat().st_mode | stat.S_IEXEC)
        print("Fixed executable permissions for mapd.")
      except Exception as e:
        print(f"Failed to set executable permissions for mapd: {e}")
        return

    try:
      print("Starting mapd...")
      subprocess.run([str(MAPD_PATH)], check=True)
    except OSError as e:
      if e.errno == 8:
        print("Exec format error encountered. Re-downloading mapd as it may be corrupted or mismatched for the architecture.")
        update_mapd()
      else:
        print(f"OSError occurred: {e}")
    except subprocess.CalledProcessError as e:
      print(f"mapd terminated with an error: {e}")
  else:
    print("mapd binary not found, waiting for network connectivity to update.")
    while not (is_url_pingable("https://github.com") or is_url_pingable("https://gitlab.com")):
      time.sleep(60)
    update_mapd()
