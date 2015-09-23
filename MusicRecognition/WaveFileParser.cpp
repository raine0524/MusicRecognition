#include <stdio.h>
#include "WaveFileParser.h"

const wave_file_t* WaveFileParser::Parse(const std::string& strFileName)
{
	FILE* fp = fopen(strFileName.c_str(), "rb");
	if (!fp)
	{
		printf("file %s open failed\n", strFileName.c_str());
		return NULL;
	}

	//read wave file header
	fread(m_stWaveFile.header.riff, sizeof(unsigned char), 4, fp);
	fread(&m_stWaveFile.header.size, sizeof(unsigned int), 1, fp);
	fread(m_stWaveFile.header.wave_flag, sizeof(unsigned char), 4, fp);
	fread(m_stWaveFile.header.fmt, sizeof(unsigned char), 4, fp);
	fread(&m_stWaveFile.header.fmt_len, sizeof(unsigned int), 1, fp);
	fread(&m_stWaveFile.header.tag, sizeof(unsigned short), 1, fp);
	fread(&m_stWaveFile.header.channels, sizeof(unsigned short), 1, fp);
	fread(&m_stWaveFile.header.samp_freq, sizeof(unsigned int), 1, fp);
	fread(&m_stWaveFile.header.byte_rate, sizeof(unsigned int), 1, fp);
	fread(&m_stWaveFile.header.block_align, sizeof(unsigned short), 1, fp);
	fread(&m_stWaveFile.header.bit_samp, sizeof(unsigned short), 1, fp);
	printf("read wave file %s header successfully\n", strFileName.c_str());

	unsigned char temp = 0;
	do { fread(&temp, sizeof(unsigned char), 1, fp); } while('d' != temp);
	m_stWaveFile.data_flag[0] = temp;
	fread(&m_stWaveFile.data_flag[1], sizeof(unsigned char), 3, fp);
	if (strncmp(m_stWaveFile.data_flag, "data", 4))
	{
		printf("Error: read flag \"data\" failed\n");
		return NULL;
	}
	fread(&m_stWaveFile.length, sizeof(unsigned int), 1, fp);
	if (m_stWaveFile.pData)
	{
		delete[] m_stWaveFile.pData;
		m_stWaveFile.pData = NULL;
	}
	m_stWaveFile.pData = new unsigned char[m_stWaveFile.length];
	fread(m_stWaveFile.pData, sizeof(unsigned char), m_stWaveFile.length, fp);
	printf("get data flag succ, and read data %d bytes\n", m_stWaveFile.length);

	fclose(fp);
	return &m_stWaveFile;
}

int WaveFileParser::Construct(const std::string& strFileName, const char* pPCMBuffer, int nBufferSize, size_t bit_per_samp, size_t samp_freq)
{
	if (strFileName.substr(strFileName.rfind(".")) != ".wav")
		return -1;

	if (!pPCMBuffer || nBufferSize <= 0 || !bit_per_samp || !samp_freq)
		return -2;

	FILE* fp = fopen(strFileName.c_str(), "wb");
	if (!fp)
	{
		printf("create file %s failed\n", strFileName.c_str());
		return -3;
	}

	unsigned int u4 = 0;
	unsigned short u2 = 0;
	fprintf(fp, "RIFF");
	u4 = 9*4+nBufferSize;
	fwrite(&u4, sizeof(unsigned int), 1, fp);		//size
	fprintf(fp, "WAVEfmt ");
	u4 = 0x00000010;
	fwrite(&u4, sizeof(unsigned int), 1, fp);
	u2 = 1;
	fwrite(&u2, sizeof(unsigned short), 1, fp);		//tag
	fwrite(&u2, sizeof(unsigned short), 1, fp);		//channels
	fwrite(&samp_freq, sizeof(size_t), 1, fp);			//sample frequency
	u4 = samp_freq*bit_per_samp/8;
	fwrite(&u4, sizeof(unsigned int), 1, fp);			//byte rate
	u2 = bit_per_samp/8;
	fwrite(&u2, sizeof(unsigned short), 1, fp);		//block align
	fwrite(&bit_per_samp, sizeof(unsigned short), 1, fp);		//bits per sample

	fprintf(fp, "data");
	fwrite(&nBufferSize, sizeof(unsigned int), 1, fp);
	fwrite(pPCMBuffer, sizeof(char), nBufferSize, fp);
	printf("construct wave file %s succ\n", strFileName.c_str());

	fclose(fp);
	return 0;
}