import os
import shutil
import datetime


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


# Path utilities
opj = os.path.join
opb = os.path.basename
opn = os.path.normpath
opd = os.path.dirname
ope = os.path.exists


#################################################################################################
#### Clean preprocessing directories
#################################################################################################
def clean(dir_prepro_raw_process, dir_prepro_fmap, dir_prepro_acpc_process, dir_prepro_orig_process,
          dir_prepro_template_process, diary_file):
    """
    Clean preprocessing directories and their contents
    """
    ct = datetime.datetime.now()
    diary = open(diary_file, "a")
    diary.write(f'\n{ct}')
    nl = '##  WARNING !!!!!!!!!!!!!!!!!!!!!       Working on step DATA CLEAN  !!!!!!!!!!!!!!!!  ## DO NOT CONFIRM THE DELETION OF DIRECTORIES IF YOU ARE NOT SURE YOU WANT TO LOSE ALL PREPROCESSED DATA!!!'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    # List of directories to remove
    list_to_remove = [dir_prepro_raw_process,
                      dir_prepro_fmap,
                      dir_prepro_acpc_process,
                      dir_prepro_orig_process,
                      dir_prepro_template_process]

    # Show what will be deleted
    print(bcolors.WARNING + "WARNING: The following directories will be PERMANENTLY DELETED:" + bcolors.ENDC)
    for directory in list_to_remove:
        if ope(directory) and os.path.isdir(directory):
            print(f"  🗑️  {directory}")
            # Show some contents to help user decide
            try:
                files = os.listdir(directory)
                if files:
                    print(f"     Contains {len(files)} files/folders")
                    if len(files) <= 5:  # Show first few items
                        for item in files[:3]:
                            print(f"       - {item}")
                        if len(files) > 3:
                            print(f"       - ... and {len(files) - 3} more")
            except:
                pass

    # Ask for confirmation
    print(bcolors.FAIL + "\n🚨 ARE YOU SURE YOU WANT TO DELETE THESE DIRECTORIES?" + bcolors.ENDC)
    print(bcolors.FAIL + "This action cannot be undone!" + bcolors.ENDC)

    response = input(bcolors.FAIL + "Type 'YES' to confirm deletion, anything else to cancel: " + bcolors.ENDC)

    if response != 'YES':
        nl = 'CANCELLED: User cancelled directory deletion.'
        print(bcolors.OKBLUE + nl + bcolors.ENDC)
        diary.write(f'\n{nl}')
        diary.close()
        return

    # Proceed with deletion
    nl = 'PROCEEDING: User confirmed directory deletion.'
    print(bcolors.OKGREEN + nl + bcolors.ENDC)
    diary.write(f'\n{nl}')

    for remove_data in list_to_remove:
        if ope(remove_data) and os.path.isdir(remove_data):
            try:
                # Use shutil.rmtree to remove directories with contents
                shutil.rmtree(remove_data)
                nl = 'SUCCESS: ' + remove_data + ' and its contents deleted!'
                print(bcolors.OKGREEN + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
            except Exception as e:
                nl = 'ERROR: Cannot delete ' + remove_data + ': ' + str(e)
                print(bcolors.FAIL + nl + bcolors.ENDC)
                diary.write(f'\n{nl}')
        else:
            nl = 'INFO: ' + remove_data + ' not found or not a directory'
            print(bcolors.WARNING + nl + bcolors.ENDC)
            diary.write(f'\n{nl}')

    # Close diary file
    diary.close()