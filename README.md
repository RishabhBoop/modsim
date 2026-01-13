# Modsim Labor (Numerische Mathematik)

> Für die deutsche Version dieser README-Datei, siehe [die deutsche Version](README-DE.md).

## Installation
This is different from the useless 'Tipps_zur_Installation'. This is using pip and venv instead of conda.
> Note: This guide has only been tested on Ubuntu 24.04 LTS (WSL2). It should work exactly the same on other Linux distributions (apart from ```apt install```). Windows users may need to adapt some commands.
### Virtual environment
Create and activate a project virtual environment (Linux/macOS):
```sh
python3 -m venv venv
source ./venv/bin/activate
```
On Windows (PowerShell):
```ps1
python -m venv venv
.\venv\Scripts\Activate.ps1
```

After activation, upgrade pip:
```sh
pip install --upgrade pip
```

### Python packages
Install the required Python packages (run these inside the activated venv):
```sh
pip install -r requirements.txt
```

### OpenModelica (Ubuntu / Debian)
1. Update packages and install prerequisites:
```sh
sudo apt-get update
sudo apt-get install -y ca-certificates curl gnupg lsb-release
```

2. Add OpenModelica GPG key:
```sh
curl -fsSL https://build.openmodelica.org/apt/openmodelica.asc | \
  sudo gpg --dearmor -o /usr/share/keyrings/openmodelica-keyring.gpg
```

3. Add the repository (uses lsb_release to determine the distro codename):
```sh
CODENAME=$(lsb_release -cs)
echo "deb [arch=amd64 signed-by=/usr/share/keyrings/openmodelica-keyring.gpg] \
  https://build.openmodelica.org/apt ${CODENAME} stable" | \
  sudo tee /etc/apt/sources.list.d/openmodelica.list
```
If lsb_release is not available, replace ${CODENAME} with your distro codename (e.g. focal, buster).

4. Install OpenModelica:
```sh
sudo apt-get update
sudo apt-get install -y openmodelica
```

## Running the code
To run the code for a specific exercise, run the corresponding Python script inside the activated virtual environment. 
Each Kapitel (K1, K2, ...) has its own folder.
All of the code for each Kapitel is contained in a single Python file (k1.py, k2.py, ...).
By default, all exercises in the script will be run.
To run only specific exercises, uncomment the corresponding function calls at the bottom of the script.
For example, to run Kapitel 4, use:
```sh
source ./venv/bin/activate  # Activate the virtual environment
python3 K4/k4.py  # Run Kapitel 4, assuming you are in the project root directory
```