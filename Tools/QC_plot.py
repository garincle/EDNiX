
import os
from nilearn import plotting
import matplotlib.pyplot as plt

opj = os.path.join
opd = os.path.dirname
ope = os.path.exists


def mosaic(img1,img2,fig):

    path_dir = opd(fig)

    if not ope(path_dir):
        os.mkdir(path_dir)
    try:
        display = plotting.plot_anat(img1,threshold='auto',display_mode='mosaic', dim=4)
        display.add_contours(img2,linewidths=.2, colors=['red'])
        display.savefig(fig)
    except:
        display = plotting.plot_anat(img2,threshold='auto',display_mode='mosaic', dim=4)
        display.savefig(fig)

    # Don't forget to close the display
    display.close()
    plt.close('all')


def mosaicB(img1, img2, fig):
    try:
        display = plotting.plot_anat(img1, display_mode='mosaic', dim=4)
        display.add_contours(img2, linewidths=.2, colors=['red'])
        display.savefig(fig)
    except:
        display = plotting.plot_anat(img2, display_mode='mosaic', dim=4)
        display.savefig(fig)

    # Don't forget to close the display
    display.close()
    plt.close('all')