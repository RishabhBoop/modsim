# Modsim Labor (Numerische Mathematik)

## Installation
Dies ist besser als die nutzlose 'Tipps_zur_Installation'. Hier wird pip und venv anstelle von conda verwendet.
> Hinweis: Diese Anleitung wurde nur unter Ubuntu 24.04 LTS (WSL2) getestet. Sie sollte auf anderen Linux-Distributionen genau gleich funktionieren (abgesehen von ```apt install```). Windows-Benutzer müssen möglicherweise einige Befehle anpassen.
### Virtuelle Umgebung
Erstellen und aktivieren Sie eine virtuelle Projektumgebung (Linux/macOS):
```sh
python3 -m venv venv
source ./venv/bin/activate
```
Unter Windows (PowerShell):
```ps1
python -m venv venv
.\venv\Scripts\Activate.ps1
```

Nach der Aktivierung, aktualisieren Sie pip:
```sh
pip install --upgrade pip
```

### Python-Pakete
Installieren Sie die erforderlichen Python-Pakete (führen Sie diese innerhalb der aktivierten venv aus):
```sh
pip install -r requirements.txt
```

### OpenModelica (Ubuntu / Debian)
1. Aktualisieren Sie die Pakete und installieren Sie die Voraussetzungen:
```sh
sudo apt-get update
sudo apt-get install -y ca-certificates curl gnupg lsb-release
```

2. Fügen Sie den OpenModelica GPG-Schlüssel hinzu:
```sh
curl -fsSL https://build.openmodelica.org/apt/openmodelica.asc | \
  sudo gpg --dearmor -o /usr/share/keyrings/openmodelica-keyring.gpg
```

3. Fügen Sie das Repository hinzu (verwendet lsb_release, um den Distro-Codenamen zu ermitteln):
```sh
CODENAME=$(lsb_release -cs)
echo "deb [arch=amd64 signed-by=/usr/share/keyrings/openmodelica-keyring.gpg] \
  https://build.openmodelica.org/apt ${CODENAME} stable" | \
  sudo tee /etc/apt/sources.list.d/openmodelica.list
```
Falls lsb_release nicht verfügbar ist, ersetzen Sie ${CODENAME} durch Ihren Distro-Codenamen (z.B. focal, buster).

4. Installieren Sie OpenModelica:
```sh
sudo apt-get update
sudo apt-get install -y openmodelica
```

## Code ausführen
Um den Code für eine bestimmte Übung auszuführen, führen Sie das entsprechende Python-Skript innerhalb der aktivierten virtuellen Umgebung aus.
Jedes Kapitel (K1, K2, ...) hat einen eigenen Ordner.
Der gesamte Code für jedes Kapitel ist in einer einzigen Python-Datei enthalten (k1.py, k2.py, ...).
Standardmäßig werden alle Übungen im Skript ausgeführt.
Um nur bestimmte Übungen auszuführen, entkommentieren Sie die entsprechenden Funktionsaufrufe am Ende des Skripts.
Um beispielsweise Kapitel 4 auszuführen, verwenden Sie:
```sh
source ./venv/bin/activate  # Aktivieren Sie die virtuelle Umgebung
python3 K4/k4.py  # Führen Sie Kapitel 4 aus, vorausgesetzt Sie befinden sich im Projektstammverzeichnis
```
