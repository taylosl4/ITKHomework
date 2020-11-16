//Template for the first part
//Averages the 21 images together
//format: ./programName outputTemplate.nii.gz

#include "iostream"
#include "string"
#include "sstream"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkAddImageFilter.h"
#include "itkDivideImageFilter.h"

const unsigned int nDims = 3 ; //only arguments are the program name and OUTPUT FILE NAME
typedef itk::Image < int, nDims > ImageType ; //this is an nii.gz

int main(int argc, const char * argv[])
{
    if( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl ;
    std::cerr << argv[0] << " dataDirectory outputGeneralTemplate" << std::endl ;
    return -1 ;
    }

    std::string names[21] = {"01","02","03","04","05","06","07","08","09","10","12","13","14","15","16","18","23","28","30","32","39"};
    std::string dir = argv[1];
    
    //Set up image reader
    typedef itk::ImageFileReader <ImageType> ImageReaderType;
    ImageReaderType::Pointer myReaderPointers[21]; //array of 21 pointers
    
    //normally we got the next part from the command line
    //in this case we will set the file name using the array of indices
    int number;
    for(int i = 0; i < 21; i++)
    {
        myReaderPointers[i]->ImageReaderType::New(); //creates new object of image reader type
        myReaderPointers[i]->SetFileName(dir + "/KKI2009-" + names[i] + "-MPRAGE.nii.gz");
        myReaderPointers[i]->Update();
    }
    
    //To average the images, we want to add them together and then divide by 21
    //Average = Sum(N) / N
    //to get the sum, we will use the AddImageFilter, which adds 2 images together
    typedef itk::AddImageFilter<ImageType,ImageType,ImageType> AddImageFilterType;
    AddImageFilterType::Pointer sumImage = AddImageFilterType::New();
    //now perform iterative addition
    sumImage->SetInput1(myReaderPointers[0]->GetOutput() );
    sumImage->SetInput2( myReaderPointers[1]->GetOutput());
    for (int i = 2; i <= 21; i++) //20 bc the last sum will be img20 + img21
    {
        sumImage->Update(); //calculates the sum
        sumImage->SetInput1(sumImage->GetOutput() );
        sumImage->SetInput2( myReaderPointers[i]->GetOutput());
    }
    sumImage->Update(); //calculates the sum
    
    //now we need to divide this output image (the output of the add filter) by 21
    typedef itk::DivideImageFilter<ImageType,ImageType,ImageType> DivideImageFilterType;
    DivideImageFilterType::Pointer avgImage = DivideImageFilterType::New();
    avgImage->SetInput(sumImage->GetOutput());
    avgImage->SetConstant(21.0);
    avgImage->Update(); //THIS IS OUR AVERAGE!!
    
    //Now to write an output
       //our output is one image that is the average of these 21 (avgImage)'
       typedef itk::ImageFileWriter<ImageType> ImageWriterType;
       ImageWriterType::Pointer myWriterPointer = ImageWriterType::New();
       myWriterPointer->SetFileName(argv[2]);
       myWriterPointer->SetInput(avgImage->GetOutput());
       myWriterPointer->Update();
    
    return 0;
}
