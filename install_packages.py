import subprocess
import sys

packages = {
    "nilearn": "0.10.4",
    "antspyx": "0.4.2",
    "nipype": "1.8.6",
    "scikit-learn": "1.5.0",  # use correct name for pip
    "numpy": "2.0.0",
    "scipy": "1.13.1",
    "nibabel": "5.2.1",
    "pandas": "2.2.2",
    "matplotlib": "3.9.2",
    "seaborn": "0.13.2",
    "pybids": "0.16.5",
    "openpyxl": "3.1.0",      # corrected version format
    "nitime": "0.11",
}

def install_with_pip(package, version):
    try:
        subprocess.check_call([
            sys.executable, "-m", "pip", "install", f"{package}=={version}"
        ])
        print(f"✅ Installed {package}=={version}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to install {package}=={version}. Error: {e}")

if __name__ == "__main__":
    for package, version in packages.items():
        install_with_pip(package, version)
