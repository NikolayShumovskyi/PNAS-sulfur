#define FILE_ALRT (1)
#define MEMORY_ALRT (2)
#define TEMP_ALRT (3)
#define VALUE_ALRT (4)

extern void options_dialog (void);
extern void too_close_dialog (int i, int j,double dist);
extern void text_error_dialog (int i);
extern void StopAlert(int i);
extern int yes(void);
extern void bond_error_dialog (int i, int j, double dist);
