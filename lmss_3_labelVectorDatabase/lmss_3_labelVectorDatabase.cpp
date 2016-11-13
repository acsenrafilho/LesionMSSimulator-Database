#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkVectorImage.h"
#include "itkComposeImageFilter.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
    if ( argc < 3 )
    {
        std::cerr << "Missing parameters. " << std::endl;
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0]
                << " labelsMainFolder sizeLesionFolder listofFilenames"
                << std::endl;
        std::cerr << "labelsMainFolder must have the output lesion size folders (5-20, 20-75,...)"
                  <<std::endl;
        std::cerr << "sizeLesionFolder informs what folder will be grouped (5-20, 20-75,75-200,200-750,750-more)"
                  <<std::endl;
        std::cerr << "listofFilenames is a list of all the file names into the sizeLesionFolder"
                  <<std::endl;
        return -1;
    }

    stringstream sizeFolder, mainFolder;
    mainFolder<<argv[1];
    sizeFolder<<argv[2];

    //    if (sizeFolder.str().c_str()!="5-20" || sizeFolder.str().c_str()!="20-75" || sizeFolder.str().c_str()!="75-200" || sizeFolder.str().c_str()!="200-750" || sizeFolder.str().c_str()!="750-more") {
    //        std::cerr<<"Wrong type of sizeFolder!"<<" ("<<sizeFolder.str().c_str()<<")"<<std::endl;
    //        return EXIT_FAILURE;
    //    }

    int numberOfInputsImages = (argc - 3 );

    const unsigned int Dimension = 3;

    typedef unsigned char                                   LabelPixelType;
    typedef itk::Image<LabelPixelType, Dimension>           ImageType;
    typedef itk::VectorImage<LabelPixelType, Dimension>     VectorImageType;

    typedef itk::ImageFileReader<ImageType> ReaderType;
    std::vector<ReaderType::Pointer> readers;

    for (int inputImg = 0; inputImg < numberOfInputsImages; ++inputImg) {
        ReaderType::Pointer reader = ReaderType::New();
        std::cout<<"Reading: "<<argv[inputImg+3]<<std::endl;
        stringstream file_path;
        file_path<<mainFolder.str().c_str()<<sizeFolder.str().c_str()<<"/"<<argv[inputImg+3];
        reader->SetFileName(file_path.str().c_str());
        reader->Update();
        readers.push_back(reader);
    }

    typedef itk::ComposeImageFilter<ImageType> ComposeType;
    typename ComposeType::Pointer compose = ComposeType::New();

    for (int n = 0; n < numberOfInputsImages; ++n) {
        compose->SetInput(n, readers[n]->GetOutput());
    }

    compose->Update();

    std::cout << compose->GetOutput()->GetNumberOfComponentsPerPixel() << std::endl;

    typedef itk::ImageFileWriter<VectorImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    stringstream output_path;
    output_path<<mainFolder.str().c_str()<<sizeFolder.str().c_str()<<"/"<<sizeFolder.str().c_str()<<"_lesions_masks.nii.gz";

    writer->SetFileName(output_path.str().c_str());
    writer->SetInput(compose->GetOutput());
    writer->Update();

    return EXIT_SUCCESS;
}

