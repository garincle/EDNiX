########### Subject loader with BIDS ##############
import os
import pandas as pd
from pathlib import Path


def Load_BIDS_to_pandas(bids_dir, modalities=None, suffixes=None, extensions=None):
    """
    Manual BIDS parser with flexible filtering options
    Completely avoids database issues and parameter limits

    Parameters
    ----------
    bids_dir : str
        Path to BIDS directory
    modalities : list, optional
        List of modalities to include: ['anat', 'func', 'dwi'] etc.
        If None, includes all standard modalities
    suffixes : list, optional
        List of file suffixes to include: ['T1w', 'T2w', 'bold', 'dwi'] etc.
        If None, includes all files
    extensions : list, optional
        List of file extensions to include: ['.nii.gz', '.nii', '.json'] etc.
        If None, includes all extensions

    Returns
    -------
    pd.DataFrame
        DataFrame with BIDS file information, ordered by subject ID and session
    """
    bids_path = Path(bids_dir)
    data = []

    # Set default modalities if not provided
    if modalities is None:
        modalities = ['anat', 'func', 'dwi', 'fmap']  # Standard BIDS modalities

    # Set default suffixes if not provided
    if suffixes is None:
        suffixes = ['T1w', 'T2w', 'bold', 'dwi', 'epi', 'phasediff', 'magnitude1', 'magnitude2']

    # Set default extensions if not provided
    if extensions is None:
        extensions = ['.nii.gz', '.nii', '.json']  # Common BIDS extensions

    print(f"Scanning BIDS directory for:")
    print(f"  Modalities: {modalities}")
    print(f"  Suffixes: {suffixes}")
    print(f"  Extensions: {extensions}")

    # Only scan valid subject directories
    subject_dirs = sorted([d for d in bids_path.glob('sub-*') if d.is_dir()])
    print(f"Found {len(subject_dirs)} subject directories")

    for subject_dir in subject_dirs:
        subject_id = subject_dir.name.replace('sub-', '')

        # Look for session directories
        session_dirs = sorted([d for d in subject_dir.glob('ses-*') if d.is_dir()])
        if not session_dirs:
            # No session structure - use subject directory directly
            session_dirs = [subject_dir]
            sessions = ['1']
        else:
            sessions = [d.name.replace('ses-', '') for d in session_dirs]

        for session_dir, session_id in zip(session_dirs, sessions):
            # Only check specified modalities
            for modality in modalities:
                modality_dir = session_dir / modality
                if modality_dir.exists():
                    # Find files with specified extensions
                    all_files = []
                    for ext in extensions:
                        # Handle both .nii.gz and simple extensions
                        if ext.startswith('.'):
                            pattern = f"*{ext}"
                        else:
                            pattern = f"*.{ext}"
                        all_files.extend(list(modality_dir.glob(pattern)))

                    for file_path in all_files:
                        file_stem = file_path.stem  # filename without extension

                        # For .nii.gz files, we need to remove .nii from stem
                        if file_path.suffix == '.gz' and file_path.stem.endswith('.nii'):
                            file_stem = file_path.stem[:-4]  # Remove .nii

                        # Extract suffix (remove sub-XXX_ses-XXX_ etc.)
                        suffix_parts = file_stem.split('_')
                        suffix = suffix_parts[-1]  # Last part is usually the suffix

                        # Get full extension
                        if file_path.suffix == '.gz' and file_path.stem.endswith('.nii'):
                            full_extension = '.nii.gz'
                        else:
                            full_extension = file_path.suffix

                        # Only include if suffix matches filter (or no filter)
                        if suffixes is None or suffix in suffixes:
                            data.append({
                                'subject': subject_id,
                                'session': session_id,
                                'suffix': suffix,
                                'extension': full_extension,
                                'datatype': modality,
                                'path': str(file_path)
                            })

    df = pd.DataFrame(data)

    # Order by subject ID and session
    if len(df) > 0:
        # Convert subject to numeric if possible for natural sorting
        try:
            # Try to convert to numeric for proper numerical ordering
            df['subject_numeric'] = pd.to_numeric(df['subject'], errors='coerce')
            # If all subjects are numeric, sort numerically
            if not df['subject_numeric'].isna().all():
                df = df.sort_values(['subject_numeric', 'session']).drop('subject_numeric', axis=1)
            else:
                # Otherwise sort alphabetically
                df = df.sort_values(['subject', 'session'])
        except:
            # Fallback: alphabetical sorting
            df = df.sort_values(['subject', 'session'])

        # Reset index for clean output
        df = df.reset_index(drop=True)

    print(f"Found {len(df)} total files matching criteria")

    # Show breakdown
    if len(df) > 0:
        print("\nFile breakdown by datatype:")
        print(df['datatype'].value_counts())
        print("\nFile breakdown by suffix:")
        print(df['suffix'].value_counts())
        print("\nFile breakdown by extension:")
        print(df['extension'].value_counts())

        # Show first few subjects in order
        print(f"\nFirst 10 subjects in order:")
        unique_subjects = df['subject'].unique()[:10]
        for subj in unique_subjects:
            subj_files = df[df['subject'] == subj]
            sessions = subj_files['session'].unique()
            print(f"  sub-{subj}: {len(subj_files)} files across sessions {sorted(sessions)}")
    else:
        print("No files found matching the specified criteria!")

    return df