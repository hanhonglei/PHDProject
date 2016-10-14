//#include "MxStdModel.h"
//
//void histoeq(MxStdModel *pMesh,float alpha /*= 0.5f*/) ;
//#include <image.h>
//#include <imagefile.h>
//
//static float round(float a) {
//  a += 0.5f;
//  return floor(a);
//}
//
//static float* histogram(const VISImage<float>& image, int nBins,
//                     float min, float max) {
//  float range = (max - min) + 1.f; 
//  float binSize = ceil(range / nBins);
// 
//  // cout << "No. of bins: " << nBins << endl;
//  // cout << "Bin size: " << binSize << endl;
//
//  int w = image.width();
//  int h = image.height();
//
//  int* bin = new int[nBins];
//  for (int i = 0; i < nBins; i++)
//    bin[i] = 0;
//
//  // float* minArr = new float[nBins];
//  // float* maxArr = new float[nBins];
//
//  /*float minMaxVal = min;
//  for (int i = 0; i < nBins; i++) {
//    minArr[i] = minMaxVal;
//    maxArr[i] = minMaxVal + binSize;
//    minMaxVal += binSize;
//  }*/
//
//  int count = 0;  
//
//  for (int i = 0; i < w; i++) {
//    for (int j = 0; j < h; j++) {
//      float pixVal = image.peek(i, j);
//
//      // correct logic, takes more iterations
//      /*for (int k = 0; k < nBins; k++) {
//        if (((pixVal >= minArr[k]) && (pixVal < maxArr[k]) && (k < (nBins - 1))) || 
//            ((pixVal >= minArr[k]) && (pixVal <= maxArr[k]) && (k == (nBins - 1)))) {
//	  bin[k] += 1;
//	  count++;
//	  break;
//	}
//      }*/
//
//      // alternate logic, takes less iterations 
//      // un-neccessary check, domain is the entire image intensity value
//      // if ((pixVal >= min) && (pixVal <= max)) {
//        int binLoc =(int)((pixVal - min)/ binSize);
//      
//        // to accomodate the max value in the range
//        if (binLoc == nBins)
//          binLoc--;
//      
//        bin[binLoc] += 1;
//	count++;
//      // }
//    }
//  }
//  
//  cout << "No. of element(s) in the bin(s): " << count << endl;
//
//  float* binFreq = new float[nBins];
//
//  // calculating the relative frequency
//  for (int i = 0; i < nBins; i++) {
//    binFreq[i] = ((float)bin[i] / (float)count);
//  }
//
//  return binFreq;  		     
//}
//
//static VISImage<float>* histoeq(const VISImage<float>& image_in, int nBins,
//                                float min, float max, float alpha, float* binFreq) {
//  float imMin = /*image_in.min()*/ min;
//  float imMax = /*image_in.max()*/ max;
//
//  float range = (imMax - imMin) + 1.f; 
//  float binSize = ceil(range / nBins);
// 
//  cout << "Bin size: " << binSize << endl;
//
//  int w = image_in.width();
//  int h = image_in.height();
//
//  VISImage<float>* image_out = new VISImage<float>(w, h); 
//  *(image_out) = 0.f;
//
//  float mnCount = 0.f;
//  for (int i = 0; i < nBins; i++) {
//    mnCount += binFreq[i];
//  }
//
//  cout << "mnCount: " << mnCount << endl;
//
//  for (int i = 0; i < w; i++) {
//    for (int j = 0; j < h; j++) {
//      float pixVal = image_in.peek(i, j);
//
//      // alternate logic, takes less iterations 
//      int binLoc =(int)((pixVal - imMin)/ binSize);
//      
//      // to accomodate the max value in the range
//      if (binLoc == nBins)
//        binLoc--;
//
//      if (binLoc > (nBins - 1)) {
//	cout << "Error in calculating binLoc\n";
//	exit(-1);
//      }	
//
//      float cdf = 0.f;
//      for (int k = 0; k <= binLoc; k++) {
//        cdf += binFreq[k];
//      }
//
//     float newPixVal = round(cdf * (max - min) + min); 
//     // float newPixVal = round(((cdf - binFreq[0]) / (mnCount - binFreq[0]))* (max - min) + min); 
//     
//     // image_out->poke(i, j) = newPixVal; 
//
//     // blending
//     image_out->poke(i, j) = alpha * newPixVal + (1.f - alpha) * pixVal;
//    }
//  }
//  
//  return image_out;  		     
//}		     
//
//int main(int argc, char **argv ) {
//  VISIm im;
//  VISImageFile imFile;
//
//  VISImage<float> image;
//  float imMin, imMax;
//  int min, max, nBins;
//  float alpha;
//
//  if ((argc > 1) && (im = imFile.read(argv[1])).isValid()) {
//    image = VISImage<float>(im);
//  }
//  else {
//    cout << "File name not entered or invalid!!!\n";
//    return -1;
//  }
//
//  imMin = image.min();
//  imMax = image.max();
//
//  if (argc == 6) {
//    nBins = atoi(argv[2]);
//    min = atoi(argv[3]);
//    max = atoi(argv[4]);
//    alpha = atof(argv[5]);
//  }
//  else if (argc == 2)  {
//    nBins = 5;
//    min = (int)imMin;
//    max = (int)imMax;
//    alpha = 1.0;
//  }
//
//  cout << "Image min and max: " << imMin << " " << imMax << endl;
//  
//  int w = image.width();
//  int h = image.height();
//
//  cout << "Image width and height: " << w << " x " << h << endl;
//
// float* binFreq =  histogram(image, nBins, /*imMin*/ min, /*imMax*/ max);
// float sum = 0.f;
//
//  for (int i = 0; i < nBins; i++) {
//    // cout << "bin[" << i << "] has frequency: " << binFreq[i] << endl; 
//    sum += binFreq[i];
//  }  
//
//  cout << "pdf's of the input image sum to: " << sum << endl;
// 
//  VISImage<float>* image_ptr =  histoeq(image, nBins, min, max, alpha, binFreq);
//  VISImage<float>& image_out = *(image_ptr);
//
//  cout << "Image min and max: " << image_out.min() << " " << image_out.max() << endl;
//  imFile.write(image_out, "!image_out.fits");
//
//  return 0;
//}
