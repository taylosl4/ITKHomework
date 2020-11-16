
import ants
dataInput = "/home/taylosl4/ANTSPyProj/KKI2009-ALL-MPRAGE"
dataOutput = "/home/taylosl4/ANTSPyProj"
from os import listdir
from os.path import isfile, join
filesList = [f for f in listdir(file_path) if isfile(join(file_path, f))]
#reading in the fixed image
fixedImage=ants.image_read("/home/taylosl4/ANTSPyProj/KKI2009-ALL-MPRAGE/KKI2009-01-MPRAGE.nii.gz");
movingImageFileList = []
movingImages = []
affineImages = []
warpedImages = []


for i in filesList: #make a list of the files that we will move through (to get a list of all of the moving images) 
    movingImageFileList.append(dataInput + "/" + i)
    movingImageList.append(f_name)
    movingImages.append((ants.image_read(dataInput + "/"+ i))) #gets filepath and filename i


for i in range(len(antsObj)): #loop through the list and run allof the registrations
    #registration: first affine with mattes MI, then symmetric normalization with MI
    transform = ants.registration(fixed=fixedImage, moving=movingImages[i],
    type_of_transform='Affine', aff_metric='mattes'))
    
    affineImages.append(ants.apply_transforms(fixed=fixedImage, moving=movingImages[i],
    transformlist=transform['fwdtransforms'] )
    
    #deformable....yes we could just run SyN but this makes it easier to seperate file writing 
    transformSyN = ants.registration(fixed=fixedImage, moving=movingImages[i],
    type_of_transform='SyNOnly', syn_metric='mattes'))
    
    warpedImages.append(ants.apply_transforms(fixed=fixedImage, moving=movingImages[i],
    transformlist=transformSyN['fwdtransforms'] ) #contains forward transforms (moving registered to fixed)
                                            
                                    
for i in range(len(file_path_m)): #writes images to files 
    ants.image_write(affineImages[i]["out"], movingImageFileList[i] + "aff.nii.gz")
    ants.image_write(warpedImages[i]["out"], movingImageFileList[i] + "warped.nii.gz")
    print('done')
