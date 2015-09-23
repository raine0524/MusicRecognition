#include <stdlib.h>
#include <string>
#include <iostream>
#include <algorithm>
#include "WaveFileParser.h"
#include "FastFourierTransform.h"

#define	PARSE_STANDARD_WAVE
//#define	PARSE_PURE_PIANO

#ifdef WIN32
	#include <Windows.h>
#else
	#include <unistd.h>
	#include <dirent.h>
	#include <sys/stat.h>
	#include <sys/types.h>
#endif

const unsigned char FastFourierTransform::BitReverseTable256[256] = 
{
	0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0, 
	0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8, 
	0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4, 
	0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC, 
	0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2, 
	0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
	0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6, 
	0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
	0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
	0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9, 
	0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
	0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
	0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3, 
	0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
	0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7, 
	0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
};

int FastFourierTransform::Radix2_FFT(const short sample[], size_t power, double amp[], int freq[], size_t K, size_t nSampFreq)
{
	if (K > (1<<power))
		return -1;

	_complex t;
	std::vector<_complex> vComplex;
	std::vector<std::pair<double, int> > vDest;
	for (size_t i = 0; i < (1<<power); i++)
	{
		t.x = sample[i];
		t.y = 0;
		vComplex.push_back(t);
	}
	ReArrangeArray(vComplex, power);
	Radix2_Dynamic(vComplex, power);
	int ret = TopKAlg(vComplex, vDest, K, nSampFreq);
	for (size_t i = 0; i < K; i++)
	{
		amp[i] = vDest[i].first;
		freq[i] = vDest[i].second;
	}
	return ret;
}

void FastFourierTransform::ReArrangeArray(std::vector<_complex>& vArray, size_t power)
{
	for (size_t i = 0; i < (1<<power); i++)
	{
		size_t c = (BitReverseTable256[i & 0xff] << 24) |
						 (BitReverseTable256[(i >> 8) & 0xff] << 16) |
						 (BitReverseTable256[(i >> 16) & 0xff] << 8) |
						 (BitReverseTable256[(i >> 24) & 0xff]);
		c >>= (32-power);
		if (c > i)
			std::swap(vArray[c], vArray[i]);
	}
}

void FastFourierTransform::Radix2_Dynamic(std::vector<_complex>& vArray, size_t power)
{
	_complex t;
	const unsigned int N = 1<<power;
	//Iteration through dyads, quadruples, octads and so on...
	for (unsigned int step = 1; step < N; step <<= 1)
	{
		//Jump to the next entry of the same transform factor
		unsigned int jump = step<<1;
		//Angle increment
		double delta = -MATH_PI/step;
		//Auxiliary sin(delta/2)
		double sine = sin(delta*0.5);

		//multiplier for trigonometric recurrence
		_complex multiplier, factor;
		multiplier.x = -2*sine*sine;
		multiplier.y = sin(delta);
		factor.x = 1.0;
		factor.y = 0.0;

		//Iteration through groups of different transform factor
		for (unsigned int group = 0; group < step; group++)
		{
			//Iteration within group
			for (unsigned int pair = group; pair < N; pair += jump)
			{
				unsigned int match = pair+step;
				_complex product;
				product.x = factor.x*vArray[match].x-factor.y*vArray[match].y;
				product.y = factor.x*vArray[match].y+factor.y*vArray[match].x;

				vArray[match].x = vArray[pair].x-product.x;
				vArray[match].y = vArray[pair].y-product.y;

				vArray[pair].x += product.x;
				vArray[pair].y += product.y;
			}
			//Successive transform factor via trigonometric recurrence
			t.x = multiplier.x*factor.x-multiplier.y*factor.y;
			t.y = multiplier.x*factor.y+multiplier.y*factor.x;
			factor.x += t.x;
			factor.y += t.y;
		}
	}
}

bool FastFourierTransform::min_heap_operator(const std::pair<double, int>& obj1, const std::pair<double, int>& obj2)
{
	return obj1.first > obj2.first;
}

int FastFourierTransform::TopKAlg(const std::vector<_complex>& vOrig, std::vector<std::pair<double, int> >& vDest, size_t K, size_t nSampFreq)
{
	const size_t N = vOrig.size();
	std::vector<std::pair<double, int> > vModule;
	for (int i = 0; i < N/2; i++)
	{
		double amp = sqrt(vOrig[i].x*vOrig[i].x+vOrig[i].y*vOrig[i].y);
		vModule.push_back(std::make_pair(amp, i*nSampFreq/N));
	}

	std::vector<std::pair<double, int> >::iterator it = vModule.begin();
	for (size_t i = 0; it != vModule.end() && i < K; it++)
	{
		if (it->second >= PIANO_FREQ_LOWER_BOUND && it->second <= PIANO_FREQ_UPPER_BOUND)
		{
			vDest.push_back(*it);
			i++;
		}
	}
	if (it == vModule.end())
		return 1;		// there're so many noisies in the frequency spectrum
	std::make_heap(vDest.begin(), vDest.end(), min_heap_operator);
	
	for (; it != vModule.end(); it++) {
		if (it->second >= PIANO_FREQ_LOWER_BOUND && it->second <= PIANO_FREQ_UPPER_BOUND && it->first > vDest.front().first) {
			std::pop_heap(vDest.begin(), vDest.end(), min_heap_operator);
			vDest[K-1] = *it;
			std::push_heap(vDest.begin(), vDest.end(), min_heap_operator);
		}
	}
	std::sort_heap(vDest.begin(), vDest.end(), min_heap_operator);
	for (std::vector<std::pair<double, int> >::iterator it = vDest.begin(); it != vDest.end(); it++)
	{
		if (it->second)
			it->first /= (N/2);
		else
			it->first /= N;
	}
	return 0;
}

FILE* g_pCalRes = NULL;

static int GetFileSize(const char* pFileName)
{
	if (!pFileName)
		return -1;

	int nFileSize;
#ifdef WIN32
	WIN32_FIND_DATA stFileInfo;
	HANDLE hFile = FindFirstFile(pFileName, &stFileInfo);
	if (INVALID_HANDLE_VALUE == hFile) {
		nFileSize = -2;
	} else {
		nFileSize =  stFileInfo.nFileSizeLow;
		FindClose(hFile);
	}
#else
	if (-1 == access(pFileName, 0)) {
		nFileSize = -3;
	} else {
		struct stat buf;
		if (-1 == stat(pFileName, &buf))
			nFileSize = -4;
		else
			nFileSize = buf.st_size;
	}
#endif
	return nFileSize;
}

void DepthFirstTraverseDir(const std::string& strRootDir, void (*pf)(void*, std::string&), void* pProcessor)
{
#ifdef WIN32
	WIN32_FIND_DATA FindData;
	HANDLE hFind = FindFirstFile(std::string(strRootDir+"*.*").c_str(), &FindData);
	if (INVALID_HANDLE_VALUE == hFind)
	{
		printf("root directory is not valid\n");
		return;
	}
	while (FindNextFile(hFind, &FindData))
	{
		if (!strcmp(FindData.cFileName, ".") || !strcmp(FindData.cFileName, ".."))
			continue;

		if (FindData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {		//directory
			DepthFirstTraverseDir(strRootDir+FindData.cFileName+"\\", pf, pProcessor);
		} else {		//file
			pf(pProcessor, strRootDir+FindData.cFileName);
			Sleep(1);
		}
	}
	FindClose(hFind);
#else
	struct dirent* ent = NULL;
	struct stat st;
	DIR* pDir = opendir(strRootDir.c_str());
	if (!pDir)
	{
		perror("open directory failed\n");
		return;
	}
	while (ent = readdir(pDir))
	{
		if (!strcmp(ent->d_name, ".") || !strcmp(ent->d_name, ".."))
			continue;

		if (-1 == stat(std::string(strRootDir+ent->d_name).c_str(), &st))
		{
			perror("stat");
			return;
		}

		if (S_ISDIR(st.st_mode)) {		//directory
			DepthFirstTraverseDir(strRootDir+ent->d_name+"/", pf, pProcessor);
		} else {		//file
			std::string strFileName = strRootDir+ent->d_name;
			pf(pProcessor, strFileName);
			usleep(1000);
		}
	}
	closedir(pDir);
#endif
}

void ProcessFile(void* pProcessor, std::string& strFileName)
{
	enum { POWER = 13, SIZE_K = 5 };
	const size_t nGroupSize = 1<<POWER;

	WaveFileParser* pWaveParser = static_cast<WaveFileParser*>(pProcessor);

#ifdef PARSE_STANDARD_WAVE
	if (strFileName.substr(strFileName.rfind(".")) != ".wav")
		return;

	const wave_file_t* wave = pWaveParser->Parse(strFileName);
	if (wave)		//parse successfully
	{
		size_t bytes = wave->header.bit_samp/8;		//number of bytes that one sampling point occupied
		/*Since the left/right channel sampling point at the same time stored alternatively,
		 *and we only need one of the specified channel sampling, so we will skip ${step} bytes
		 *whereby channels=1 means that there exist just one channel, =2 means have two channels.
		 */
		size_t step = bytes*(wave->header.channels-1), nGroups = wave->length/nGroupSize/bytes/wave->header.channels;
		unsigned char* pData = wave->pData;

		short sample[nGroupSize];
		double amp[SIZE_K];
		int freq[SIZE_K];
		short nSamplePoint;		//temporary variable

		for (size_t i = 0; i < nGroups; i++)
		{
			for (size_t j = 0; j < nGroupSize; j++, pData += step)		//construct a set of sample
			{
				nSamplePoint = 0;
				for (size_t k = 0; k < bytes; k++, pData++)		// construct a sample point
					nSamplePoint |= (*pData)<<(8*k);
				sample[j] = nSamplePoint;
			}

			//run FastFourierTransform on this group of sample
			int nRetCode = FastFourierTransform::Radix2_FFT(sample, POWER, amp, freq, SIZE_K, wave->header.samp_freq);
			if (!nRetCode)
			{
				fprintf(g_pCalRes, "%s: ",strFileName.c_str());
				for (int index = 0; index < SIZE_K; index++)
					fprintf(g_pCalRes, "%lf,%d\t\t", amp[index], freq[index]);
				fprintf(g_pCalRes, "\n\n");
			}
		}
	}
#endif

#ifdef PARSE_PURE_PIANO
	if (strFileName.substr(strFileName.rfind(".")) != ".txt")
		return;

	int nFileSize = GetFileSize(strFileName.c_str());
	if (nFileSize <= 0)
		return;

	char* pPCMBuffer = new char[nFileSize];
	FILE* pFile = fopen(strFileName.c_str(), "rb");
	fread(pPCMBuffer, 1, nFileSize, pFile);
	fclose(pFile);

	size_t bit_per_sample = 16, samp_freq = /*16000*/44100;
	size_t nGroups = nFileSize/nGroupSize/(bit_per_sample/8);
	short* pData = (short*)pPCMBuffer;

	std::vector<double> vSample;
	std::vector<std::pair<double, int> > vResult;
	short nSamplePoint;

	for (size_t i = 0; i < nGroups; i++)
	{
		for (size_t j = 0; j < nGroupSize; j++, pData++)
		{
			nSamplePoint = *pData;
			vSample.push_back(nSamplePoint);
		}

		//run FastFourierTransform on this group of sample
		FastFourierTransform::Radix2_FFT(vSample, POWER, vResult, SIZE_K, samp_freq);
		fprintf(g_pCalRes, "%s: ",strFileName.c_str());
		for (int index = 0; index < SIZE_K; index++)
			fprintf(g_pCalRes, "%lf,%d\t\t", vResult[index].first, vResult[index].second);
		fprintf(g_pCalRes, "\n\n");
		vSample.clear();
	}

	strFileName.replace(strFileName.rfind("."), strlen(".txt"), ".wav");
	pWaveParser->Construct(strFileName, pPCMBuffer, nFileSize, bit_per_sample, samp_freq);
	delete[] pPCMBuffer;
#endif
}

//test case
double physical_signal(double t)
{
	return 2+3*cos(2*MATH_PI*50*t-MATH_PI*30/180)+1.5*cos(2*MATH_PI*75*t+MATH_PI*90/180)+5*sin(2*MATH_PI*440*t+MATH_PI*60/180)+7.5*cos(2*MATH_PI*365*t+MATH_PI*45/180);
}

int main(int argc, char* argv[])
{
#if 1
	enum { POWER = 12, SIZE_K = 5 };
	short test[1<<POWER];
	double amp[SIZE_K];
	int freq[SIZE_K];
	double nTimeUnit = 1.0/(1<<POWER);
	for (size_t i = 0; i < (1<<POWER); i++)
		test[i] = physical_signal(i*nTimeUnit);
	
	int nRetCode = FastFourierTransform::Radix2_FFT(test, POWER, amp, freq, SIZE_K, 1<<POWER);
	if (!nRetCode)
	{
		for (size_t i = 0; i < SIZE_K; i++)
			std::cout<<amp[i]<<'\t'<<freq[i]<<std::endl;
	}
#else
	g_pCalRes = fopen("CalculateResult.txt", "w");
#ifdef PARSE_STANDARD_WAVE
	std::string strDir = "E:\\res\\MusicRecognition\\piano_test\\";
#endif
#ifdef PARSE_PURE_PIANO
	std::string strDir = "..\\piano_test\\";
#endif
	WaveFileParser WaveParser;
	DepthFirstTraverseDir(strDir, ProcessFile, &WaveParser);
	fclose(g_pCalRes);
#endif

	system("pause");
	return EXIT_SUCCESS;
}