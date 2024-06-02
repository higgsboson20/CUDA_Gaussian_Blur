
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846

/**This function is to create the gaussian kernel matrix 
float** create_kernel_matrix(size_t kernel_order, uint32_t sigma) {
    float kernel_matrix[kernel_order][kernel_order];
    float normalizer = 1 / (2 * 3.14 * sigma * sigma);
    for (size_t x = 0; x < kernel_order; x++) {
        for (size_t y = 0; y < kernel_order; y++) {
            kernel_matrix[x][y] = normalizer * exp(-1 * (x * x + y * y) / (2 * sigma * sigma));
        }
    }
    return kernel_matrix;
}
**/

int main(int argc, char * argv[])
{   
    /**To ensure there are three arguments: <input_pgm>, <output_pgm> and <sigma>*/
    if (argc != 4) {
         fprintf(stderr, "Usage: ./gaussian_blur_serial <input_pgm> <output_pgm> <sigma>\n");
    }

    /**Parse the command line to recieve the input file, output file name and the sigma value*/
    const char* input_file = argv[1]; //input binary PGM file
    const char* output_file = argv[2]; //output binary PGM file
    uint32_t sigma = atoi(argv[3]); //sigma value

    char temp[200]; // temporary string to hold a copy of the input filename
    uint32_t width = 0;

    // Copy the filename to a temporary string
    strcpy(temp, input_file);

    // Find the position of the last underscore
    char *underscore_pos = strrchr(temp, '_');
    if (underscore_pos != NULL) {
        // Find the position of the dot
        char *dot_pos = strchr(underscore_pos, '.');
        if (dot_pos != NULL) {
            // Null-terminate the string at the dot position to isolate the number
            *dot_pos = '\0';
            // Convert the substring to an integer
            width = atoi(underscore_pos + 1);
        }
    }
    printf("%d\n", width);

    /**Open the input binary PGM file*/
    FILE* file = fopen(input_file, "rb");
    if (file == NULL) {
        perror("Cannot open file");
        exit(1);
    }

    /**Validate the bounds of sigma*/
    if (sigma <= 0) {
        fprintf(stderr, "Sigma value must be greater than 0\n");
        return EXIT_FAILURE;
    }
    
    /**Create the gaussian kernel matrix*/
    uint32_t kernel_order = sigma * 6; //order of the kernel matrix
    uint32_t half_order = kernel_order / 2;

    if (kernel_order % 2 == 0) kernel_order += 1; //make it odd to account for edges

    //float** kernel_matrix = create_kernel_matrix(kernel_order, sigma);
    float kernel_matrix[kernel_order][kernel_order];
    float normalizer = 1 / (2 * PI * sigma * sigma);
    float sum = 0.0;

    //fills in the gaussian kernel matrix
    for (size_t x = 0; x < kernel_order; x++) {
        for (size_t y = 0; y < kernel_order; y++) {
            uint32_t x_offset = x - half_order;
            uint32_t y_offset = y - half_order;
            kernel_matrix[x][y] = normalizer * exp(-1.0 * (x_offset * x_offset + y_offset * y_offset) / (2.0 * sigma * sigma));
            sum += kernel_matrix[x][y];
        }
    }

    //normalizes the gaussian kernel matrix to the total sum
    for (size_t i = 0; i < kernel_order; i++) {
        for (size_t j = 0; j < kernel_order; j++) {
            kernel_matrix[i][j] /= sum;
            printf("%.8f\t", kernel_matrix[i][j]);
        }
        printf("\n");
    }

    return 0;

} 