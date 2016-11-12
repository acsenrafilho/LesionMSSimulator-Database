#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
    if ( argc < 4 )
    {
        std::cerr << "Missing parameters. " << std::endl;
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0]
                << " imagesFolder imagePrefix minimumLesionSize(pixels) outputLesionFolder"
                << std::endl;
        std::cerr<<"The imagesFolder must contain the images files: <imagePrefix>_T1.nii.gz, <imagePrefix>_T2.nii.gz, <imagePrefix>_FLAIR.nii.gz, <imagePrefix>_FA.nii.gz, <imagePrefix>_ADC.nii.gz, <imagePrefix>_label.nii.gz"
                << std::endl;
        std::cerr << "outputLesionFolder must have the output lesion size folders (0-20, 20-75,...)"
                  <<std::endl;
        return -1;
    }
    stringstream images_folder, imagePrefix, output_folder;
    images_folder<<argv[1];
    imagePrefix<<argv[2];
    output_folder<<argv[4];
    int lesionMinSize = atoi(argv[3]);

    const unsigned int Dimension = 3;

    typedef float                                    PixelType;
    typedef unsigned char                            LabelPixelType;
    typedef itk::Image<PixelType, Dimension>         ImageType;
    typedef itk::Image<LabelPixelType, Dimension>    LabelImageType;

    typedef itk::ImageFileReader<ImageType> ReaderType;
    typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
    typedef itk::ImageFileWriter<LabelImageType> WriterType;

    ReaderType::Pointer imageT1 = ReaderType::New();
    ReaderType::Pointer imageT2 = ReaderType::New();
    ReaderType::Pointer imageFLAIR = ReaderType::New();
    ReaderType::Pointer imageFA = ReaderType::New();
    ReaderType::Pointer imageADC = ReaderType::New();
    LabelReaderType::Pointer label = LabelReaderType::New();

    stringstream t1_path,t2_path,flair_path,fa_path, adc_path, label_path;
    t1_path<<images_folder.str().c_str()<<imagePrefix.str().c_str()<<"_T1.nii.gz";
    t2_path<<images_folder.str().c_str()<<imagePrefix.str().c_str()<<"_T2.nii.gz";
    flair_path<<images_folder.str().c_str()<<imagePrefix.str().c_str()<<"_FLAIR.nii.gz";
    fa_path<<images_folder.str().c_str()<<imagePrefix.str().c_str()<<"_FA.nii.gz";
    adc_path<<images_folder.str().c_str()<<imagePrefix.str().c_str()<<"_ADC.nii.gz";
    label_path<<images_folder.str().c_str()<<imagePrefix.str().c_str()<<"_label.nii.gz";

    imageT1->SetFileName(t1_path.str().c_str());
    imageT2->SetFileName(t2_path.str().c_str());
    imageFLAIR->SetFileName(flair_path.str().c_str());
    imageFA->SetFileName(fa_path.str().c_str());
    imageADC->SetFileName(adc_path.str().c_str());
    label->SetFileName(label_path.str().c_str());

    try
    {
//        imageT1->Update();
//        imageT2->Update();
//        imageFLAIR->Update();
//        imageFA->Update();
//        imageADC->Update();
        label->Update();
    }
    catch ( itk::ExceptionObject &err)
    {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return -1;
    }

    //Label treatment - Removing and organizing the lesion database
    //Spliting labels
    typedef itk::ConnectedComponentImageFilter<LabelImageType, LabelImageType> ConnectedLabelType;
    typename ConnectedLabelType::Pointer connLabel = ConnectedLabelType::New();
    connLabel->SetInput(label->GetOutput());
    connLabel->Update();

    cout<<"Number of connected lesions: "<<connLabel->GetObjectCount()<<endl;

    //Cleaning areas lower than 3 mm3 of volume
    typedef itk::RelabelComponentImageFilter<LabelImageType, LabelImageType>      RelabelerType;
    typename RelabelerType::Pointer relabel = RelabelerType::New();
    relabel->SetInput(connLabel->GetOutput());
    relabel->SetSortByObjectSize(true);
    relabel->SetMinimumObjectSize(lesionMinSize);
    relabel->Update();

    int nLesion = relabel->GetNumberOfObjects();
    cout<<"Number of considered lesions: "<<nLesion<<endl;
    cout<<"Size of the lesions: ";
    for (unsigned int size = 0; size < nLesion; ++size) {
        if (size<nLesion-1) {
            cout<<relabel->GetSizeOfObjectInPixels(size)<<",";
        }else{
            cout<<relabel->GetSizeOfObjectInPixels(size)<<endl;
        }
    }

    typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType>   BinaryThresholdType;
    typedef itk::ImageFileWriter<LabelImageType>                      WriterType;
    typename BinaryThresholdType::Pointer binary = BinaryThresholdType::New();
    typename WriterType::Pointer writer = WriterType::New();

    binary->SetInput(relabel->GetOutput());

    for (int l = 0; l <= nLesion; ++l) {
        binary->SetInsideValue(1);
        binary->SetOutsideValue(0);
        //Thresholding each lesion from the original volume
        binary->SetLowerThreshold(l);
        binary->SetUpperThreshold(l+1);
        binary->Update();

        if (relabel->GetSizeOfObjectInPixels(l)>=750) {
            //Saving the lesion map: 750-more mm3
            stringstream output_file;
            output_file<<output_folder.str().c_str()<<"750-more";
            if (l<10) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-00"<<l;
            }else if (l>=10 & l<100) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-0"<<l;
            }else if (l>=100 & l<1000) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-"<<l;
            }

            cout<<"Saving label - "<<l<<endl;
            output_file<<".nii.gz";
            writer->SetInput(binary->GetOutput());
            writer->SetFileName(output_file.str().c_str());
            writer->Update();
        }else if (relabel->GetSizeOfObjectInPixels(l)<750 && relabel->GetSizeOfObjectInPixels(l)>=200) {
            //Saving the lesion map: 200-750 mm3
            stringstream output_file;
            output_file<<output_folder.str().c_str()<<"200-750";
            if (l<10) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-00"<<l;
            }else if (l>=10 & l<100) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-0"<<l;
            }else if (l>=100 & l<1000) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-"<<l;
            }

            cout<<"Saving label - "<<l<<endl;
            output_file<<".nii.gz";
            writer->SetInput(binary->GetOutput());
            writer->SetFileName(output_file.str().c_str());
            writer->Update();
        }else if (relabel->GetSizeOfObjectInPixels(l)<200 && relabel->GetSizeOfObjectInPixels(l)>=75) {
            //Saving the lesion map: 75-200 mm3
            stringstream output_file;
            output_file<<output_folder.str().c_str()<<"75-200";
            if (l<10) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-00"<<l;
            }else if (l>=10 & l<100) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-0"<<l;
            }else if (l>=100 & l<1000) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-"<<l;
            }

            cout<<"Saving label - "<<l<<endl;
            output_file<<".nii.gz";
            writer->SetInput(binary->GetOutput());
            writer->SetFileName(output_file.str().c_str());
            writer->Update();
        }else if (relabel->GetSizeOfObjectInPixels(l)<75 && relabel->GetSizeOfObjectInPixels(l)>=20) {
            //Saving the lesion map: 20-75 mm3
            stringstream output_file;
            output_file<<output_folder.str().c_str()<<"20-75";
            if (l<10) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-00"<<l;
            }else if (l>=10 & l<100) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-0"<<l;
            }else if (l>=100 & l<1000) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-"<<l;
            }

            cout<<"Saving label - "<<l<<endl;
            output_file<<".nii.gz";
            writer->SetInput(binary->GetOutput());
            writer->SetFileName(output_file.str().c_str());
            writer->Update();
        }else if (relabel->GetSizeOfObjectInPixels(l)<20 && relabel->GetSizeOfObjectInPixels(l)>=5) {
            //Saving the lesion map: 5-20 mm3
            stringstream output_file;
            output_file<<output_folder.str().c_str()<<"5-20";
            if (l<10) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-00"<<l;
            }else if (l>=10 & l<100) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-0"<<l;
            }else if (l>=100 & l<1000) {
                output_file<<"/"<<imagePrefix.str().c_str()<<"-lesion-"<<l;
            }

            cout<<"Saving label - "<<l<<endl;
            output_file<<".nii.gz";
            writer->SetInput(binary->GetOutput());
            writer->SetFileName(output_file.str().c_str());
            writer->Update();
        }



    }



    //Extract information from the images
    typedef itk::ImageFileWriter<ImageType>         WriterImageType;
    typename WriterImageType::Pointer writerImage = WriterImageType::New();

    //T1
    typename ImageType::Pointer contrastT1 = ImageType::New();
    contrastT1->CopyInformation(imageT1->GetOutput());
    contrastT1->Allocate();




    return EXIT_SUCCESS;
}

