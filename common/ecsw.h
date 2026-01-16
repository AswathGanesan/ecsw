enum ecsw_mode {
    ECSW_TEST,
    ECSW_BENCHMARK
};

extern void setup(const enum ecsw_mode mode);
extern void teardown(void);
extern void output(char *message);
extern unsigned long long  performance_count(void);
extern char* performance_count_unit(void);
