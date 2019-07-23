#define ASSERT_LEVEL 1 

inline void _assert(int assert_level, bool assertion, const char* exp, const char* txt, 
					const char* file, int line, ...)
{
	
  	#if ASSERT_LEVEL > 0

		// a. check assertion
		if(assert_level >= ASSERT_LEVEL || assertion == true){
			return;
		}

	    // b. print general info
	    FILE* log_file = fopen("log_assert.txt","w");
	    fprintf(log_file, "%s failed \nfile %s, line %d\n", exp, file, line);

	    // c. print user specified info
	    fprintf(log_file, "messsage:\n");
	    fprintf(log_file, "  ");
		va_list args;
		va_start (args, line);  
		vfprintf (log_file, txt, args);
	    va_end (args);

	    // d. clean up
	    fclose(log_file);
	    mexErrMsgTxt("assertion not true, see log_assert.txt");

  	#endif

}
#define assert(LEVEL,EXP,TXT,...) _assert(LEVEL,EXP, #EXP, TXT, __FILE__, __LINE__, ##__VA_ARGS__)