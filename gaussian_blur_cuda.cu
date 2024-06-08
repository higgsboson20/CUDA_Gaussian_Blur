#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846

#define DIV_ROUND_UP(n, d)  (((n) + (d) - 1) / (d))

#define cuda_check(ret) _cuda_check((ret), __FILE__, __LINE__)

inline void _cuda_check(cudaError_t ret, const char *file, int line)
{
        if (ret != cudaSuccess) {
                fprintf(stderr, "CudaErr: %s (%s:%d)\n", cudaGetErrorString(ret), file, line);
                exit(1);
        }
}

__global__ void convolution_kernel(int kernel_order, float * kernel_matrix, unsigned char * image_in, unsigned char * image_out, size_t width, size_t height)
{
    size_t row = blockIdx.y*blockDim.y + threadIdx.y;
    size_t col = blockIdx.x*blockDim.x + threadIdx.x;

    int dc = (kernel_order - 1)/2;

    float conv = 0, sum = 0;

    for(size_t k_i = 0; k_i < kernel_order; k_i++){
        for(size_t k_j = 0; k_j < kernel_order; k_j++){
            int ii = row + k_i - dc;
            int jj = col + k_j - dc;

            // Handle boundaries by clamping to the nearest edge pixel
            if (ii < 0) ii = 0;
            if (ii >= height) ii = height - 1;
            if (jj < 0) jj = 0;
            if (jj >= width) jj = width - 1;

            conv += image_in[width*ii + jj]*kernel_matrix[kernel_order*k_i + k_j];
            sum += kernel_matrix[kernel_order*k_i + k_j];

        }
    }

    if(sum > 0)
        conv /= sum;

    image_out[row*width + col] = (unsigned char)conv;

}

int main(int argc, char * argv[])
{
    /**To ensure there are three arguments: <input_pgm>, <output_pgm> and <sigma>*/
    if (argc != 4) {
         fprintf(stderr, "Usage: ./gaussian_blur_serial <input_pgm> <output_pgm> <sigma>\n");
    }


    // =================================================
    // Fill initial image buffer, initialize out image buffer, and create kernel matrix
    // =================================================

    /**Parse the command line to recieve the input file, output file name and the sigma value*/
    const char* input_file = argv[1]; //input binary PGM file
    const char* output_file = argv[2]; //output binary PGM file
    float sigma = atof(argv[3]); //sigma value

    //char temp[200]; // temporary string to hold a copy of the input filename
    size_t width, height;
    unsigned char max_value;

    /**Open the input binary PGM file*/
    FILE* file = fopen(input_file, "rb");
    if (file == NULL) {
        perror("Cannot open file");
        exit(1);
    }

    if (fscanf(file, "P5\n%ld %ld\n%hhd\n", &width, &height, &max_value) != 3) {
        fprintf(stderr, "Error reading image header\n");
        fclose(file);
        return EXIT_FAILURE;
    }

    unsigned char * image = (unsigned char *)malloc(width*height); // image binary data
    unsigned char * blurred_image = (unsigned char *)malloc(width*height); // blurred image
    /** Push image into buffer */
    uint32_t seek_dist = 3 + (int)log10(height) + (int)log10(width) + 4 + 4; // |firstLine| + |secondLine| + |thirdLine|
    fseek(file, seek_dist, SEEK_SET);

    if (fread(image, 1, height * width, file) != height * width) {
        fprintf(stderr, "Error reading image data\n");
        free(image);
        fclose(file);
        return 1;
    }

    /**Validate the bounds of sigma*/
    if (sigma <= 0) {
        fprintf(stderr, "Sigma value must be greater than 0\n");
        exit(1);
    }

    /**Create the gaussian kernel matrix*/
    uint32_t kernel_order = ceil(sigma * 6); //order of the kernel matrix
    uint32_t half_order = kernel_order / 2;

    if (kernel_order % 2 == 0) kernel_order += 1; //make it odd to account for edges

    // float** kernel_matrix = create_kernel_matrix(kernel_order, sigma);
    float * kernel_matrix  = (float *)malloc(sizeof(float) * kernel_order * kernel_order);
    float normalizer = 1 / (2 * PI * sigma * sigma);
    float sum = 0.0;

    //fills in the gaussian kernel matrix
    for (size_t x = 0; x < kernel_order; x++) {
        for (size_t y = 0; y < kernel_order; y++) {
            uint32_t x_offset = x - half_order;
            uint32_t y_offset = y - half_order;
            kernel_matrix[kernel_order*x + y] = normalizer * exp(-1.0 * (x_offset * x_offset + y_offset * y_offset) / (2.0 * sigma * sigma));
            sum += kernel_matrix[kernel_order*x + y];
        }
    }

    /** normalizes the gaussian kernel matrix to the total sum */
    for (size_t i = 0; i < kernel_order; i++) {
        for (size_t j = 0; j < kernel_order; j++) {
            kernel_matrix[kernel_order*i + j] /= sum;
        }
    }

    // =================================================
    //
    // =================================================


    // =================================================
    // initialize CUDA business
    // =================================================

    unsigned char * image_in, * image_out;
    float * kernel_matrix_d;

    cuda_check(cudaMalloc(&kernel_matrix_d, kernel_order*kernel_order*sizeof(float)));
    cuda_check(cudaMalloc(&image_in, width*height));
    cuda_check(cudaMalloc(&image_out, width*height));

    cuda_check(cudaMemcpy(kernel_matrix_d, kernel_matrix, kernel_order*kernel_order*sizeof(float), cudaMemcpyHostToDevice));
    cuda_check(cudaMemcpy(image_in, image, width*height, cudaMemcpyHostToDevice));

    dim3 block_dim(32,32);
    dim3 grid_dim(DIV_ROUND_UP(width, block_dim.x), DIV_ROUND_UP(height, block_dim.y));
    convolution_kernel<<<grid_dim, block_dim>>>(kernel_order, kernel_matrix_d, image_in, image_out, width, height);

    cuda_check(cudaPeekAtLastError());
    cuda_check(cudaDeviceSynchronize());

    cuda_check(cudaMemcpy(blurred_image, image_out, width*height, cudaMemcpyDeviceToHost));

    cuda_check(cudaFree(image_in));
    cuda_check(cudaFree(image_out));


    // =================================================
    //
    // =================================================

    /** output file */
    FILE * out = fopen(output_file, "wb");
    if(!out){
        fprintf(stderr, "Error with opening file!\n");
        exit(1);
    }
    fprintf(out, "P5\n%ld %ld\n%d\n", width, height, max_value);

    fwrite(blurred_image, sizeof(char), width*height, out);

    fclose(out);
    fclose(file);

    free(image);
    free(kernel_matrix);
    free(blurred_image);



    return 0;

}
