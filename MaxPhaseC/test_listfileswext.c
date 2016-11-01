#include "maxphase.h"
#include <stdio.h>

int main(){
	
	char *extArr[] = {"a","b","abc"};
	
	char *dirName = "testListFileswExtDir";
		
	size_t nFiles, maxFileNameLen, lpc;
	
	char **fileList;
	
	fileList = listfileswext(extArr[0], dirName, &nFiles, &maxFileNameLen);
    printf("Found %zu files with extension .%s\n",nFiles,extArr[0]);
	printf("Length of longest file name: %zu\n",maxFileNameLen);
	for(lpc = 0; lpc < nFiles; lpc++){
		printf("%s\n", fileList[lpc]);
	    /* Free up dynamically allocated memory */
		free(fileList[lpc]);
	}
	/* Free up dynamically allocated memory */
	free(fileList);
	
	fileList = listfileswext(extArr[1], dirName, &nFiles, &maxFileNameLen);
    printf("Found %zu files with extension .%s\n",nFiles,extArr[1]);
	printf("Length of longest file name: %zu\n",maxFileNameLen);
	for(lpc = 0; lpc < nFiles; lpc++){
		printf("%s\n", fileList[lpc]);
	    /* Free up dynamically allocated memory */
		free(fileList[lpc]);
	}
	/* Free up dynamically allocated memory */
	free(fileList);
	
	fileList = listfileswext(extArr[2], dirName, &nFiles, &maxFileNameLen);
    printf("Found %zu files with extension .%s\n",nFiles,extArr[2]);
	printf("Length of longest file name: %zu\n",maxFileNameLen);
	for(lpc = 0; lpc < nFiles; lpc++){
		printf("%s\n", fileList[lpc]);
	    /* Free up dynamically allocated memory */
		free(fileList[lpc]);
	}
	/* Free up dynamically allocated memory */
	free(fileList);
	
}