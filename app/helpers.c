#define CFG_FILE_DELIMITER "= "

// Prints statistics about app execution
void printStats(
        const struct timeval appExecutionStarted, 
        const struct timeval appExecutionEnded, 
        const int currentProcess, 
        const int totalProcesses) {

    const long difference = (appExecutionEnded.tv_sec - appExecutionStarted.tv_sec) * 1000.0f 
        + (appExecutionEnded.tv_usec - appExecutionStarted.tv_usec) / 1000.0f;

    printf("\nProcess %d/%d has completed in %ldms\n", currentProcess + 1, totalProcesses, difference);
}

// Returns maximum complex number
long double max(long double *array, const int size) {
    long double lastMaximum = array[0];

    for (int i = 1 ; i < size ; i++) {
        if (array[i] > lastMaximum) {
            lastMaximum = array[i];
        }
    }

    return lastMaximum;
}

// Prints a vector of complex numbers
void printListFileForSage(const long double complex *vector, const int size, const int iterationId, const int processId) {
        char processIdString[sizeof(int)];
        sprintf(processIdString, "%d", processId);
        
        char fullPath[1024];
        strcpy(fullPath, "result-");
        strcat(fullPath, processIdString);
        strcat(fullPath, ".txt");

        FILE *file = fopen(fullPath, "a");
        
        fprintf(file, "\n\nlist%d = [", iterationId);
        
        for (int i = 0 ; i <= size ; i++) {
            fprintf(file, "%.8Lf", cabsl(vector[i]));

            if (i != size) {
                fprintf(file, ", ");
            }
        }

		fprintf(file, "]");
        
        fclose(file);
}

// Read configuration by process id
struct configuration getConfigurationBy(const int processId) {
    struct configuration cfg;

    char processIdString[sizeof(int)];
    sprintf(processIdString, "%d", processId);


    char fullPath[1024];
    strcpy(fullPath, "conf/");
    strcat(fullPath, processIdString);
    strcat(fullPath, ".cfg");

    FILE *file = fopen(fullPath, "r");

    char line[1024];
    int lineId = 0;

    while (fgets(line, sizeof(line), file) != NULL) {

        char *configValue;
        configValue = strstr(line, CFG_FILE_DELIMITER);
        configValue = configValue + strlen(CFG_FILE_DELIMITER);

        switch (lineId) {
            case 0:
                cfg.N = atoi(configValue);
                break;

            case 1:
                cfg.Tau = atof(configValue);
                break;

            case 2:
                cfg.T = atof(configValue);
                break;

            case 3:
                cfg.alpha = atof(configValue);
                break;

            case 4:
                cfg.delta = atof(configValue);
                break;
        }

        lineId++;
    }

    fclose(file);

    return cfg;
}