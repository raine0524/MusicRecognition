#ifndef	FASTFOURIERTRANSFORM_H_
#define	FASTFOURIERTRANSFORM_H_
#include <vector>
#include <utility>
#include <math.h>

#define	MATH_PI		3.14159265358979323846
#define	PIANO_FREQ_LOWER_BOUND	25
#define	PIANO_FREQ_UPPER_BOUND	4250

class FastFourierTransform
{
public:
	/* \brief fast fourier transform master procedure
	 * @sample[]: sampling set with complex form, the real part of each sampling point indicate
	 * the current time point while the imaginary part indicate the amplitude of the wave
	 * @array.size() == (1<<power), the size of array had better to be 2 of integer power
	 * @amp[] & freq[]: the amplitudes and frequencies of the K sine waves with the maximum amplitude,
	 * which were the OUTPUT arrays as the FFT's result
	 * @K: indicates the size of OUTPUT array
	 * @nSampleFreq: sample frequency
	 */
	static int Radix2_FFT(const short sample[], size_t power, double amp[], int freq[], size_t K, size_t nSampFreq);

private:
	/* \brief get the K maximum frequencies in the frequency domain
	 * @vOrig: complex array after processed by the Radix2_FFT
	 * @vDest: contains K sine waves, it->first indicates the amplitude while it->second indicates the frequency
	 * @K: K pairs
	 * @nSampleFreq: sample frequency
	 */
	static int TopKAlg(const std::vector<_complex>& vOrig, std::vector<std::pair<double, int> >& vDest, size_t K, size_t nSampFreq);

	/* \brief a pre-process operator used in the Radix-2 FFT, processed object is
	 * \complex vector, and in this application the size of array is restricted to 
	 * \2 of integer power as well as less than 2^16=65536, but the function
	 * \still have the general ability in 32-bit processor.
	 */
	static void ReArrangeArray(std::vector<_complex>& vArray, size_t power);

	/* \brief use the butterfly operator process the scramble array
	 */
	static void Radix2_Dynamic(std::vector<_complex>& vArray, size_t power);

	/* \brief callback function used in the std::make_heap to make minimum heap
	 */
	static bool min_heap_operator(const std::pair<double, int>& obj1, const std::pair<double, int>& obj2);

private:
	static const unsigned char BitReverseTable256[256];
};

#endif		//FASTFOURIERTRANSFORM_H_