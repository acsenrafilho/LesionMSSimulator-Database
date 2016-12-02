#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkMaskNegatedImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"

#include <iostream>
#include <fstream>

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
        std::cerr<<"The imagesFolder must contain the images files: <imagePrefix>_T1.nii.gz, <imagePrefix>_T2.nii.gz, <imagePrefix>_PD.nii.gz, <imagePrefix>_FLAIR.nii.gz, <imagePrefix>_FA.nii.gz, <imagePrefix>_ADC.nii.gz, <imagePrefix>_label.nii.gz"
                << std::endl;
        std::cerr << "outputLesionFolder must have the output lesion size folders (50-100, 100-500,...)"
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
    typedef unsigned int                            LabelPixelType;
    typedef itk::Image<PixelType, Dimension>         ImageType;
    typedef itk::Image<LabelPixelType, Dimension>    LabelImageType;

    typedef itk::ImageFileReader<ImageType> ReaderType;
    typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
    typedef itk::ImageFileWriter<LabelImageType> WriterType;

    ReaderType::Pointer imageT1 = ReaderType::New();
    ReaderType::Pointer imageT2 = ReaderType::New();
    ReaderType::Pointer imageFLAIR = ReaderType::New();
    ReaderType::Pointer imagePD = ReaderType::New();
    ReaderType::Pointer imageFA = ReaderType::New();
    ReaderType::Pointer imageADC = ReaderType::New();
    LabelReaderType::Pointer label = LabelReaderType::New();

    stringstream t1_path,t2_path,flair_path,pd_path, fa_path, adc_path, label_path;
    t1_path<<images_folder.str().c_str()<<imagePrefix.str().c_str()<<"_T1.nii.gz";
    t2_path<<images_folder.str().c_str()<<imagePrefix.str().c_str()<<"_T2.nii.gz";
    flair_path<<images_folder.str().c_str()<<imagePrefix.str().c_str()<<"_FLAIR.nii.gz";
    pd_path<<images_folder.str().c_str()<<imagePrefix.str().c_str()<<"_PD.nii.gz";
    fa_path<<images_folder.str().c_str()<<imagePrefix.str().c_str()<<"_FA.nii.gz";
    adc_path<<images_folder.str().c_str()<<imagePrefix.str().c_str()<<"_ADC.nii.gz";
    label_path<<images_folder.str().c_str()<<imagePrefix.str().c_str()<<"_label.nii.gz";

    imageT1->SetFileName(t1_path.str().c_str());
    imageT2->SetFileName(t2_path.str().c_str());
    imageFLAIR->SetFileName(flair_path.str().c_str());
    imagePD->SetFileName(pd_path.str().c_str());
    imageFA->SetFileName(fa_path.str().c_str());
    imageADC->SetFileName(adc_path.str().c_str());
    label->SetFileName(label_path.str().c_str());

    try
    {
        imageT1->Update();
        imageT2->Update();
        imageFLAIR->Update();
        imagePD->Update();
        imageFA->Update();
        imageADC->Update();
        label->Update();
    }
    catch ( itk::ExceptionObject &err)
    {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return -1;
    }

    //Rescale images to 8-bits
    //T1
    typedef itk::RescaleIntensityImageFilter<ImageType,ImageType>       RescalerType;
    typename RescalerType::Pointer imageT18b = RescalerType::New();
    imageT18b->SetInput(imageT1->GetOutput());
    imageT18b->SetOutputMinimum(0.0);
    imageT18b->SetOutputMaximum(255.0);
    imageT18b->Update();

    //T2
    typedef itk::RescaleIntensityImageFilter<ImageType,ImageType>       RescalerType;
    typename RescalerType::Pointer imageT28b = RescalerType::New();
    imageT28b->SetInput(imageT2->GetOutput());
    imageT28b->SetOutputMinimum(0.0);
    imageT28b->SetOutputMaximum(255.0);
    imageT28b->Update();

    //T2-FLAIR
    typedef itk::RescaleIntensityImageFilter<ImageType,ImageType>       RescalerType;
    typename RescalerType::Pointer imageFLAIR8b = RescalerType::New();
    imageFLAIR8b->SetInput(imageFLAIR->GetOutput());
    imageFLAIR8b->SetOutputMinimum(0.0);
    imageFLAIR8b->SetOutputMaximum(255.0);
    imageFLAIR8b->Update();

    //PD
    typedef itk::RescaleIntensityImageFilter<ImageType,ImageType>       RescalerType;
    typename RescalerType::Pointer imagePD8b = RescalerType::New();
    imagePD8b->SetInput(imagePD->GetOutput());
    imagePD8b->SetOutputMinimum(0.0);
    imagePD8b->SetOutputMaximum(255.0);
    imagePD8b->Update();

    //FA
    typedef itk::RescaleIntensityImageFilter<ImageType,ImageType>       RescalerType;
    typename RescalerType::Pointer imageFA8b = RescalerType::New();
    imageFA8b->SetInput(imageFA->GetOutput());
    imageFA8b->SetOutputMinimum(0.0);
    imageFA8b->SetOutputMaximum(255.0);
    imageFA8b->Update();

    //ADC
    typedef itk::RescaleIntensityImageFilter<ImageType,ImageType>       RescalerType;
    typename RescalerType::Pointer imageADC8b = RescalerType::New();
    imageADC8b->SetInput(imageADC->GetOutput());
    imageADC8b->SetOutputMinimum(0.0);
    imageADC8b->SetOutputMaximum(255.0);
    imageADC8b->Update();

    //******************************************************************
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

        if (relabel->GetSizeOfObjectInPixels(l)>=5000) {
            //Saving the lesion map: 5000-more mm3
            stringstream output_file;
            output_file<<output_folder.str().c_str()<<"5000-more";
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
        }else if (relabel->GetSizeOfObjectInPixels(l)<5000 && relabel->GetSizeOfObjectInPixels(l)>=1000) {
            //Saving the lesion map: 1000-5000 mm3
            stringstream output_file;
            output_file<<output_folder.str().c_str()<<"1000-5000";
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
        }else if (relabel->GetSizeOfObjectInPixels(l)<1000 && relabel->GetSizeOfObjectInPixels(l)>=500) {
            //Saving the lesion map: 500-1000 mm3
            stringstream output_file;
            output_file<<output_folder.str().c_str()<<"500-1000";
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
        }else if (relabel->GetSizeOfObjectInPixels(l)<500 && relabel->GetSizeOfObjectInPixels(l)>=100) {
            //Saving the lesion map: 100-500 mm3
            stringstream output_file;
            output_file<<output_folder.str().c_str()<<"100-500";
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
        }else if (relabel->GetSizeOfObjectInPixels(l)<100 && relabel->GetSizeOfObjectInPixels(l)>=50) {
            //Saving the lesion map: 50-100 mm3
            stringstream output_file;
            output_file<<output_folder.str().c_str()<<"50-100";
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

    //******************************************************************
    //Extract information from the images
    fstream measureDatabase;
    stringstream measure_filepath;
    measure_filepath<<output_folder.str().c_str()<<"/"<<imagePrefix.str().c_str()<<"_measuresDatabase.csv";
    cout<<"Saving table in: "<<measure_filepath.str().c_str()<<endl;
    measureDatabase.open(measure_filepath.str().c_str(), ios::out);
    measureDatabase<<",T1,,,,,T2,,,,,T2-FLAIR,,,,,PD,,,,,FA,,,,,ADC,,,,\n";
    measureDatabase<<"pat_label,contrast,mean,std,min,max,contrast,mean,std,min,max,contrast,mean,std,min,max,contrast,mean,std,min,max,contrast,mean,std,min,max,contrast,mean,std,min,max\n";

    //White matter mask (contrast calculation)
    LabelReaderType::Pointer wm_mask = LabelReaderType::New();
    stringstream wm_mask_path;
    wm_mask_path<<output_folder.str().c_str()<<"/MNI152_T1_1mm_brain_wm.nii.gz";
    wm_mask->SetFileName(wm_mask_path.str().c_str());
    wm_mask->Update();

    //Masking lesion
    typedef itk::MaskNegatedImageFilter<LabelImageType, LabelImageType>     MaskNegateType;
    typename MaskNegateType::Pointer maskN = MaskNegateType::New();
    maskN->SetInput(wm_mask->GetOutput());

    //Statistic
    typedef itk::LabelStatisticsImageFilter<ImageType, LabelImageType>      StatisticType;
    typename StatisticType::Pointer stat = StatisticType::New();

    for (int l = 1; l <= nLesion; ++l) {
        binary->SetInsideValue(1);
        binary->SetOutsideValue(0);
        //Thresholding each lesion from the original volume
        binary->SetLowerThreshold(l);
        binary->SetUpperThreshold(l+1);
        binary->Update();

        maskN->SetMaskImage(binary->GetOutput());
        maskN->Update();

        cout<<"Extracting information from label: "<<l<<endl;
        double contrast=0.0;
        if (relabel->GetSizeOfObjectInPixels(l)>=5) {
            //Extracting information
            //T1
            //Measure: Contrast
            stat->SetInput(imageT18b->GetOutput());
            stat->SetLabelInput(binary->GetOutput());
            stat->Update();
            contrast=stat->GetMean(1);

            stat->SetInput(imageT18b->GetOutput());
            stat->SetLabelInput(maskN->GetOutput());
            stat->Update();
            contrast=contrast/stat->GetMean(3);

            if (l<10) {
                measureDatabase<<imagePrefix.str().c_str()<<"-lesion-00"<<l<<",";
            }else if (l>=10 & l<100) {
                measureDatabase<<imagePrefix.str().c_str()<<"-lesion-0"<<l<<",";
            }else if (l>=100 & l<1000) {
                measureDatabase<<imagePrefix.str().c_str()<<"-lesion-"<<l<<",";
            }

            measureDatabase<<contrast<<",";


            //Measure: Mean,Sigma,Max and Min (from lesion)
            stat->SetInput(imageT18b->GetOutput());
            stat->SetLabelInput(binary->GetOutput());
            stat->Update();

            measureDatabase<<stat->GetMean(1)<<",";
            measureDatabase<<stat->GetSigma(1)<<",";
            measureDatabase<<stat->GetMinimum(1)<<",";
            measureDatabase<<stat->GetMaximum(1)<<",";

            //*****
            //T2
            //*****
            //Measure: Contrast
            stat->SetInput(imageT28b->GetOutput());
            stat->SetLabelInput(binary->GetOutput());
            stat->Update();
            contrast=stat->GetMean(1);

            stat->SetInput(imageT28b->GetOutput());
            stat->SetLabelInput(maskN->GetOutput());
            stat->Update();
            contrast=contrast/stat->GetMean(3);

            measureDatabase<<contrast<<",";

            //Measure: Mean,Sigma,Max and Min (from lesion)
            stat->SetInput(imageT28b->GetOutput());
            stat->SetLabelInput(binary->GetOutput());
            stat->Update();

            measureDatabase<<stat->GetMean(1)<<",";
            measureDatabase<<stat->GetSigma(1)<<",";
            measureDatabase<<stat->GetMinimum(1)<<",";
            measureDatabase<<stat->GetMaximum(1)<<",";

            //*****
            //T2-FLAIR
            //*****
            //Measure: Contrast
            stat->SetInput(imageFLAIR8b->GetOutput());
            stat->SetLabelInput(binary->GetOutput());
            stat->Update();
            contrast=stat->GetMean(1);

            stat->SetInput(imageFLAIR8b->GetOutput());
            stat->SetLabelInput(maskN->GetOutput());
            stat->Update();
            contrast=contrast/stat->GetMean(3);

            measureDatabase<<contrast<<",";

            //Measure: Mean,Sigma,Max and Min (from lesion)
            stat->SetInput(imageFLAIR8b->GetOutput());
            stat->SetLabelInput(binary->GetOutput());
            stat->Update();

            measureDatabase<<stat->GetMean(1)<<",";
            measureDatabase<<stat->GetSigma(1)<<",";
            measureDatabase<<stat->GetMinimum(1)<<",";
            measureDatabase<<stat->GetMaximum(1)<<",";

            //*****
            //PD
            //*****
            //Measure: Contrast
            stat->SetInput(imagePD8b->GetOutput());
            stat->SetLabelInput(binary->GetOutput());
            stat->Update();
            contrast=stat->GetMean(1);

            stat->SetInput(imagePD8b->GetOutput());
            stat->SetLabelInput(maskN->GetOutput());
            stat->Update();
            contrast=contrast/stat->GetMean(3);

            measureDatabase<<contrast<<",";

            //Measure: Mean,Sigma,Max and Min (from lesion)
            stat->SetInput(imagePD8b->GetOutput());
            stat->SetLabelInput(binary->GetOutput());
            stat->Update();

            measureDatabase<<stat->GetMean(1)<<",";
            measureDatabase<<stat->GetSigma(1)<<",";
            measureDatabase<<stat->GetMinimum(1)<<",";
            measureDatabase<<stat->GetMaximum(1)<<",";

            //*****
            //FA
            //*****
            //Measure: Contrast
            stat->SetInput(imageFA8b->GetOutput());
            stat->SetLabelInput(binary->GetOutput());
            stat->Update();
            contrast=stat->GetMean(1);

            stat->SetInput(imageFA8b->GetOutput());
            stat->SetLabelInput(maskN->GetOutput());
            stat->Update();
            contrast=contrast/stat->GetMean(3);

            measureDatabase<<contrast<<",";

            //Measure: Mean,Sigma,Max and Min (from lesion)
            stat->SetInput(imageFA8b->GetOutput());
            stat->SetLabelInput(binary->GetOutput());
            stat->Update();

            measureDatabase<<stat->GetMean(1)<<",";
            measureDatabase<<stat->GetSigma(1)<<",";
            measureDatabase<<stat->GetMinimum(1)<<",";
            measureDatabase<<stat->GetMaximum(1)<<",";

            //*****
            //ADC
            //*****
            //Measure: Contrast
            stat->SetInput(imageADC8b->GetOutput());
            stat->SetLabelInput(binary->GetOutput());
            stat->Update();
            contrast=stat->GetMean(1);

            stat->SetInput(imageADC8b->GetOutput());
            stat->SetLabelInput(maskN->GetOutput());
            stat->Update();
            contrast=contrast/stat->GetMean(3);

            measureDatabase<<contrast<<",";

            //Measure: Mean,Sigma,Max and Min (from lesion)
            stat->SetInput(imageADC8b->GetOutput());
            stat->SetLabelInput(binary->GetOutput());
            stat->Update();

            measureDatabase<<stat->GetMean(1)<<",";
            measureDatabase<<stat->GetSigma(1)<<",";
            measureDatabase<<stat->GetMinimum(1)<<",";
            measureDatabase<<stat->GetMaximum(1)<<",";

            //End of line
            measureDatabase<<"\n";
        }
    }
    measureDatabase.close();
    return EXIT_SUCCESS;
}

