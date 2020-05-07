// For benchmarking purposes only
static uint64_t flip_count = 0;

// Flip buffering variables
static uint32_t flip_word = 0;
static int flip_pos = 0;
int flip(void)
{
	if(flip_pos == 0)
	{
		flip_word = genrand_int32();
		flip_pos = 32;
	}
	flip_count++;
	flip_pos--;
	return (flip_word & (1 << flip_pos)) >> flip_pos;
}

inline uint32_t algFDR(unsigned int n)
{
	uint32_t v = 1, c = 0;
	while(true)
	{
		v = v << 1;
		c = (c << 1) + flip();
		if(v >= n)
		{
			if(c < n) return c;
			else
			{
				v = v - n;
				c = c - n;
			}
		}
	}
}
