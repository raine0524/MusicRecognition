#ifndef	WAVEFILEPARSER_H_
#define	WAVEFILEPARSER_H_
#include <string.h>
#include <string>

//wave header
typedef struct tagwave_header_t
{
	unsigned char riff[4];					//资源交换文件标志
	unsigned int size;							//从下个地址开始到文件结尾的字节数
	unsigned char wave_flag[4];		//wave文件标识
	unsigned char fmt[4];					//波形格式标识
	unsigned int fmt_len;					//过滤字节（一般为0000 0010H）
	unsigned short tag;						//格式种类，值为1时，标识PCM线性编码
	unsigned short channels;			//通道数，单声道为1，双声道为2
	unsigned int samp_freq;				//采样频率
	unsigned int byte_rate;				//数据传输率（每秒字节=采样频率*每个样本字节数）
	unsigned short block_align;		//块对齐字节数=channels*bit_samp/8
	unsigned short bit_samp;			//bits per sample（又称量化位数）

	tagwave_header_t()
	{
		memset(this, 0, sizeof(tagwave_header_t));
	}
} wave_header_t;

typedef struct tagwave_file_t
{
	wave_header_t header;			//header
	char data_flag[4];						//数据标识符
	unsigned int length;				//采样数据总数
	unsigned char* pData;			//data

	tagwave_file_t()
	{
		memset(this, 0, sizeof(tagwave_file_t));
	}
	~tagwave_file_t()
	{
		if (pData)
		{
			delete[] pData;
			pData = NULL;
		}
	}
} wave_file_t;

class WaveFileParser
{
public:
	WaveFileParser() {}
	~WaveFileParser() {}

	const wave_file_t* Parse(const std::string& strFileName);
	int Construct(const std::string& strFileName, const char* pPCMBuffer, int nBufferSize, size_t bit_per_samp, size_t samp_freq);

private:
	wave_file_t m_stWaveFile;
};

#endif		//WAVEFILEPARSER_H_