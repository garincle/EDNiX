import os
import glob
import json
import subprocess
from os.path import join as opj
import ants

# FONCTION SPÉCIFIQUE POUR MOTION CORRECTION + TRANSFORMATIONS
def apply_motion_correction_and_transforms(input_4d, reference_image, output_4d,
                                           mat_files_pattern, additional_transformations,
                                           TR, sing_afni, temp_dir=None, interpolator='linear',
                                           metadata_sources=None, description=None):
    """
    Applique motion correction (matrices individuelles) + autres transformations
    Version Python ANTs (pas de commandes shell)
    """

    # Trouver toutes les matrices de motion correction
    mat_files = sorted(glob.glob(mat_files_pattern))
    if not mat_files:
        raise ValueError(f"Aucun fichier .mat trouvé avec le pattern: {mat_files_pattern}")

    print(f"Found {len(mat_files)} motion correction matrices")

    # Créer le répertoire temporaire
    if temp_dir is None:
        temp_dir = opj(os.path.dirname(output_4d), 'tmp')
    os.makedirs(temp_dir, exist_ok=True)

    commands_executed = []

    # 1. Split de l'image 4D avec AFNI
    split_prefix = opj(temp_dir, 'dummy.nii.gz')
    split_cmd = f'{sing_afni}3dTsplit4D -overwrite -prefix {split_prefix} {input_4d}'

    print("Splitting 4D volume with 3dTsplit4D...")
    subprocess.run(split_cmd, shell=True, check=True)
    commands_executed.append(split_cmd)

    # 2. Trouver tous les fichiers splités
    WW = sorted(glob.glob(opj(temp_dir, 'dummy*.nii.gz')))

    if len(WW) != len(mat_files):
        print(f"⚠️ Attention: {len(WW)} volumes mais {len(mat_files)} matrices de motion")

    # 3. Charger l'image de référence
    ref_img = ants.image_read(reference_image)

    # 4. Appliquer les transformations à chaque volume avec ANTs Python
    transformed_files = []

    for j, (input_file, mat_file) in enumerate(zip(WW, mat_files)):
        # Générer le nom de fichier de sortie
        if j < 10:
            output_temp = opj(temp_dir, f'warp1.00{j}.nii.gz')
        elif j < 100:
            output_temp = opj(temp_dir, f'warp1.0{j}.nii.gz')
        else:
            output_temp = opj(temp_dir, f'warp1.{j}.nii.gz')

        # Charger l'image d'entrée
        moving_img = ants.image_read(input_file)

        # Construire la liste des transformations pour CE volume
        transform_list = [mat_file] + additional_transformations

        # Appliquer avec ANTs Python (pas de commande shell)
        transformed_img = ants.apply_transforms(
            fixed=ref_img,
            moving=moving_img,
            transformlist=transform_list,
            interpolator=interpolator
        )

        # Sauvegarder l'image transformée
        ants.image_write(transformed_img, output_temp)

        # Enregistrer la commande "virtuelle"
        transform_cmd = f"ants.apply_transforms with {len(transform_list)} transformations"
        commands_executed.append(transform_cmd)

        print(f"Transformed volume {j + 1}/{len(WW)} with {os.path.basename(mat_file)}")
        transformed_files.append(output_temp)

    # 5. Recréer le 4D avec AFNI
    files_string = ' '.join(transformed_files)
    tcat_cmd = f'{sing_afni}3dTcat {files_string} -tr {TR} -overwrite -prefix {output_4d}.nii.gz'

    print("Creating 4D volume with 3dTcat...")
    subprocess.run(tcat_cmd, shell=True, check=True)
    commands_executed.append(tcat_cmd)

    # 6. Métadonnées
    if metadata_sources is None:
        metadata_sources = [input_4d, reference_image] + mat_files + additional_transformations

    metadata = {
        "Sources": metadata_sources,
        "Description": description or f"Motion correction + {len(additional_transformations)} transformations",
        "Command": commands_executed,
        "Parameters": {
            "interpolator": interpolator,
            "TR": TR,
            "motion_matrices": len(mat_files),
            "additional_transformations": len(additional_transformations)
        }
    }

    json_file = output_4d.replace('.nii.gz', '.json')
    with open(json_file, "w") as outfile:
        json.dump(metadata, outfile, indent=3)

    print(f"✓ Motion correction + transformations terminé: {output_4d}")
    return output_4d