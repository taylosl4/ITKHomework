//Template for the affine registration
//finds affine transformation matrix of 21 images
//performs registration, averages
//format: ./programName outputTemplate.nii.gz

#include "iostream"
#include "string"
#include "sstream"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkAddImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkAffineTransform.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkCommand.h"

//setup
const unsigned int nDims = 3 ;
typedef itk::Image < int, nDims > ImageType ;

class OptimizerIterationCallback : public itk::Command
{
public:
  // itk set up such that things like pointers and new() are available for this class
  //OptimizerIterationCallback::Pointer myCallback = OptimizerIterationCallback::New() ;
  // standard things we need in every itk-inheriting class
  typedef OptimizerIterationCallback Self ;
  typedef itk::Command Superclass ;
  typedef itk::SmartPointer<OptimizerIterationCallback> Pointer ;
  itkNewMacro(OptimizerIterationCallback);

  // specific typedefs we need for our observer
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType ;
  typedef const OptimizerType * OptimizerPointerType ;

  // if i want to change things in my caller
 void
 Execute(itk::Object * caller, const itk::EventObject & event)
  {
    Execute((const itk::Object *)caller, event);
  }
 
  // if i am just observing the caller (no changes)
  void
  Execute(const itk::Object * caller, const itk::EventObject & event)
  {
    // somehow get my hands on the optimizer
    // caller is the optimizer, but it's of the wrong type, but i know it is an optimizer so i can cast it to the optimizer type
    OptimizerPointerType optimizer = dynamic_cast < OptimizerPointerType > ( caller )  ;
    std::cout << optimizer->GetCurrentIteration() << " " << optimizer->GetValue() << std::endl ;
  }
};

class RegistrationIterationCallback : public itk::Command
{
public:
  // itk set up such that things like pointers and new() are available for this class
  //OptimizerIterationCallback::Pointer myCallback = OptimizerIterationCallback::New() ;
  // standard things we need in every itk-inheriting class
  typedef RegistrationIterationCallback Self ;
  typedef itk::Command Superclass ;
  typedef itk::SmartPointer<RegistrationIterationCallback> Pointer ;
  itkNewMacro(RegistrationIterationCallback);

  // specific typedefs we need for our observer
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType ;
  typedef OptimizerType * OptimizerPointerType ;
  
  typedef itk::MultiResolutionImageRegistrationMethod < ImageType, ImageType > RegistrationMethodType ;
  typedef RegistrationMethodType * RegistrationPointerType ;

  // if i want to change things in my caller
 void
 Execute(itk::Object * caller, const itk::EventObject & event)
  {
    // somehow get my hands on the registration method
    // caller is the registration method
    RegistrationPointerType registration = dynamic_cast < RegistrationPointerType > ( caller ) ;
    std::cout << "Level: " << registration->GetCurrentLevel () << std::endl ;

    OptimizerPointerType optimizer = dynamic_cast < OptimizerPointerType > ( registration->GetModifiableOptimizer() ) ;

    optimizer->SetMaximumStepLength ( optimizer->GetMaximumStepLength() * 0.5 ) ;
    /*
    if ( registration->GetCurrentLevel() == 0 )
      optimizer->SetMaximumStepLength ( 0.25 ) ;
    else if ( registration->GetCurrentLevel() == 1 )
      optimizer->SetMaximumStepLength ( 0.125 ) ;
    else
      optimizer->SetMaximumStepLength ( 0.125 ) ;
    */
  }
 
  // if i am just observing the caller (no changes)
  void
  Execute(const itk::Object * caller, const itk::EventObject & event)
  {
    // nothing
  }
};

int main ( int argc, char * argv[] )
{
  // Verify command line arguments
  if( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl ;
    std::cerr << argv[0] << "outputRegisteredImageFile" << std::endl ;
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
    
  // Registration!

  // set up typedefs (don't forget to include the header files)
  //  typedef itk::ImageRegistrationMethod < ImageType, ImageType > RegistrationMethodType ;
  typedef itk::MultiResolutionImageRegistrationMethod < ImageType, ImageType > RegistrationMethodType ;
  typedef itk::AffineTransform < double, 3 > AffineTransformType ; // leaving at default values, could skip
  typedef itk::MeanSquaresImageToImageMetric < ImageType, ImageType > MetricType ;
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType ; // no template arguments
  typedef itk::LinearInterpolateImageFunction < ImageType > InterpolatorType ;
  typedef itk::ResampleImageFilter < ImageType, ImageType > ResampleFilterType ;

  // declare the variables
  RegistrationMethodType::Pointer registrationMethod = RegistrationMethodType::New() ;
  AffineTransformType::Pointer transform = AffineTransformType::New ();
  MetricType::Pointer metric = MetricType::New () ;
  OptimizerType::Pointer optimizer = OptimizerType::New() ;
  InterpolatorType::Pointer interpolator = InterpolatorType::New() ;

    typedef itk::ImageFileWriter < ImageType > ImageWriterType ;
    ImageWriterType::Pointer resampledWriter[21];
    ImageReaderType::Pointer resampledReader[21]; //array of 21 pointers

  // connect the pipeline
    for (int i = 1; i< 21; i++)
    {
        registrationMethod->SetMovingImage ( myReaderPointers[i]->GetOutput() ) ;
        registrationMethod->SetFixedImage ( myReaderPointers[0]->GetOutput() ) ;
        registrationMethod->SetOptimizer ( optimizer ) ;
        registrationMethod->SetMetric ( metric ) ;
        registrationMethod->SetInterpolator ( interpolator ) ;
        registrationMethod->SetTransform ( transform ) ;
        
      // set up the relevant parameters
      optimizer->MinimizeOn () ;
      optimizer->SetNumberOfIterations ( 20 ) ;

      std::cout << "Min: " << optimizer->GetMinimumStepLength () << std::endl ;
      std::cout << "Max: " << optimizer->GetMaximumStepLength () << std::endl ;
      std::cout << "Current: " << optimizer->GetCurrentStepLength () << std::endl ;

      optimizer->SetMinimumStepLength ( 0 ) ;
      optimizer->SetMaximumStepLength ( 0.125 ) ; // TODO: might want to adjust this some more
      transform->SetIdentity () ;
      registrationMethod->SetInitialTransformParameters ( transform->GetParameters() ) ;
     
      // set up the callback function
      OptimizerIterationCallback::Pointer myCallback = OptimizerIterationCallback::New();
      optimizer->AddObserver(itk::IterationEvent(), myCallback);

      // set up the second callback for the registration method
      RegistrationIterationCallback::Pointer myRegCallback = RegistrationIterationCallback::New() ;
      registrationMethod->AddObserver ( itk::IterationEvent(), myRegCallback ) ;

      registrationMethod->SetNumberOfLevels ( 3 ) ;
      registrationMethod->SetFixedImageRegion ( myReaderPointers[0]->GetOutput()->GetLargestPossibleRegion() ) ;

      // run the registration
      // TODO: put this in a try-catch block
      registrationMethod->Update() ;

      // why did it stop?
      std::cout << optimizer->GetStopConditionDescription () << std::endl ;

      // apply the transform we get from the registration
      ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New() ;
      resampleFilter->SetInput ( myReaderPointers[i]->GetOutput() ) ;
      // transform that we are applying
      resampleFilter->SetTransform ( transform ) ;
      // set the grid to where the fixed image is, just in case it moved too far
      resampleFilter->SetReferenceImage ( myReaderPointers[0]->GetOutput() ) ;
      resampleFilter->UseReferenceImageOn () ;

      // update
      resampleFilter->Update() ;
        //first write resampled images to file
      resampledWriter[i]->ImageWriterType::New();
      resampledWriter[i]->SetFileName(dir + "temp" + names[i] + "nii.gz" );
      resampledWriter[i]->SetInput(resampleFilter->GetOutput());
      resampledWriter[i]->Update();
        
      
      //now READ this resampled imaged
      resampledReader[i]->ImageReaderType::New();
      resampledReader[i]->SetFileName(dir + "temp" + names[i] + "nii.gz");
      resampledReader[i]->Update();
        
    }
    
    //now we need to average these images
    //lets read them from the reader
    typedef itk::AddImageFilter<ImageType,ImageType,ImageType> AddImageFilterType;
    AddImageFilterType::Pointer sumImage = AddImageFilterType::New();
    //now perform iterative addition
    sumImage->SetInput1(resampledReader[0]->GetOutput() );
    sumImage->SetInput2( resampledReader[1]->GetOutput());
    for (int i = 2; i <= 21; i++) //20 bc the last sum will be img20 + img21
    {
        sumImage->Update(); //calculates the sum
        sumImage->SetInput1(sumImage->GetOutput() );
        sumImage->SetInput2( resampledReader[i]->GetOutput());
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
   ImageWriterType::Pointer myWriterPointer = ImageWriterType::New();
   myWriterPointer->SetFileName(argv[2]);
   myWriterPointer->SetInput(avgImage->GetOutput());
   myWriterPointer->Update();

  // Done.
  return 0 ;
}

