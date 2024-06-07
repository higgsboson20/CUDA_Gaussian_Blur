
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#define PI 3.14159265358979323846

int main(int argc, char * argv[])
{   
    /**To ensure there are three arguments: <input_pgm>, <output_pgm> and <sigma>*/
    if (argc != 4) {
         fprintf(stderr, "Usage: ./gaussian_blur_serial <input_pgm> <output_pgm> <sigma>\n");
    }

    /**Parse the command line to recieve the input file, output file name and the sigma value*/
    const char* input_file = argv[1]; //input binary PGM file
    const char* output_file = argv[2]; //output binary PGM file
    float sigma = atof(argv[3]); //sigma value

    /**Validate the bounds of sigma and kernel*/
    if (sigma <= 0) {
        fprintf(stderr, "Sigma value must be greater than 0\n");
        return EXIT_FAILURE;
    }

    uint32_t width, height, max_value = 0;

    /**Open the input binary PGM file*/
    FILE* file = fopen(input_file, "rb");
    if (!file) {
        perror("Cannot open file");
        exit(1);
    }

    // Read PGM header
    if (fscanf(file, "P5\n%d %d\n%d\n", &width, &height, &max_value) != 3) {
        fprintf(stderr, "Error reading image header\n");
        fclose(file);
        return 1;
    }

    printf("%d\n", width);
    printf("%d\n", height);
    printf("%d\n", max_value);

    //get the image matrix
    unsigned char* input_image = malloc(height * width);
    unsigned char* blurred_image =  malloc(height * width);
    uint32_t seek_dist = 3 + 2*(int)log10(width) + 4 + 4; // |firstLine| + |secondLine| + |thirdLine|
    fseek(file, seek_dist, SEEK_SET);
    if (fread(input_image, sizeof(unsigned char), height * width, file) != height * width) {
        fprintf(stderr, "Error reading image data\n");
        free(input_image);
        fclose(file);
        return 1;
    }
/** 
     for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++) {
            printf("%.8f\t", input_matrix[i][j]);
        }
        printf("\n");
    }
    **/
    printf("height: %d\n", height);
    printf("width: %d\n", width);

    /**Create the gaussian kernel matrix*/
    uint32_t kernel_order = sigma * 6; //order of the kernel matrix
    
    /**Validate the bounds of the kernel matrix*/
    if (kernel_order > width || kernel_order > height) {
        fprintf(stderr, "The kernel matrix should not be bigger than the input image size \n");
        return EXIT_FAILURE;
    }

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

    /**Perform the convolution of the kernel matrix onto the input image*/
    for (size_t i = 0; i < height; i++) {
        for (size_t j = 0; j < width; j++) {
            float pixel_value = 0.0;

            for (size_t k = 0; k < kernel_order; k++) {
                for (size_t l = 0; l < kernel_order; l++) {
                    size_t x = j - half_order + l;
                    size_t y = i - half_order + k;  

                    if (x >= 0 && x < width && y >= 0 && y < height) {
                        pixel_value += input_image[width * y + x] * kernel_matrix[k][l];
                    }
                }
            }
            blurred_image[i * width + j]= (unsigned char) pixel_value;
        }
    }

    /**Write to the PGM file */
    FILE* out = fopen(output_file, "wb");
    if (!out) {
        perror("Cannot open file");
        exit(1);
    }

    //write PGM header
    fprintf(out, "P5\n%d %d\n%d\n", width, height, max_value);
    fwrite(blurred_image, sizeof(unsigned char), height * width, out);
    fclose(out);
    free(input_image);
    free(blurred_image);

    return 0;
} 