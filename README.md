# Modsim Labor (Numerische Mathematik)

## Installation
Tnis is different from the useless 'Tipps_zur_Installation'. This is using pip and venv instead of conda.
> Note: This guide has only been tested on Ubuntu 24.04 LTS (WSL2). It should work exactly the same on other Linux distributions. Windows users may need to adapt some commands.
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