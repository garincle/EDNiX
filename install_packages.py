import subprocess

# Dictionary of packages and their required versions
packages = {
    "shutil": "built-in (no version)",
    "subprocess": "built-in (no version)",
    "os": "built-in (no version)",
    "glob": "built-in (no version)",
    "json": "built-in (no version)",
    "nilearn": "0.10.4",
    "ants": "0.4.2",
    "nipype": "1.8.6",
    "sklearn": "1.5.0",
    "copy": "built-in (no version)",
    "numpy": "2.0.0",
    "scipy": "1.13.1",
    "nibabel": "5.2.1",
    "math": "built-in (no version)",
    "pandas": "2.2.2",
    "matplotlib": "3.9.2",
    "pingouin": "0.5.4",
    "seaborn": "0.13.2",
    "bids": "0.16.5",
    "openpyxl": "3.1.00",
    "torch": "2.4.1",
    "nitime": "0.11",
    "numba": "0.61.0rc1"}

def install_with_conda(package, version):
    try:
        subprocess.check_call(["conda", "install", "-y", f"{package}={version}"])
        print(f"Successfully installed {package}={version} with conda")
    except subprocess.CalledProcessError:
        print(f"Failed to install {package}={version} with conda. Trying with pip...")
        return False
    return True

def install_with_pip(package, version):
    try:
        subprocess.check_call(["pip", "install", f"{package}=={version}"])
        print(f"Successfully installed {package}=={version} with pip")
    except subprocess.CalledProcessError as e:
        print(f"Failed to install {package}=={version} with pip. Error: {e}")

if __name__ == "__main__":
    for package, version in packages.items():
        if "built-in" not in version and "not installed" not in version:
            if not install_with_conda(package, version):
                install_with_pip(package, version)
