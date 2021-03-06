#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N  22  		// количество базисных переменных 
#define M  4   		// количество переменных состояния 
#define BRANCHES  5 // количество ветвей нормировнного дерева
#define HORDES  4	// количество хорд

// параметры элементов схемы
#define C13 1e-7 				// c1
#define C14 1e-7 				// c2

#define L 1e-3
#define E 10.
#define R 1e3

// VD11
#define C_d11  2e-12			// c3
#define R_d11  1e6				// r1
#define It_d11  1e-12
#define MFt_d11  0.026
#define Rb_d11  20.				// r2

template <class T>
void init_array(T* array, int size_x, int size_y){
	int i, j;	

	for(j = 0; j < size_y; j++)
		for(i = 0 ; i < size_x ; i++){
			array[size_x * j + i] = 0.;
		}
}

template <class T>
void print_array(T* array, int size_x, int size_y){
	int i, j;	
  
	for(j = 0; j < size_y; j++){
		for(i = 0 ; i < size_x ; i++){
			printf("%5.1lf ", array[size_x * j + i]);
		}
		printf("\n\n");
	}
}

template <class T>
void print_array_to_file(T* array, int size_x, int size_y, FILE* fp){
	int i, j;	
  
	for(j = 0; j < size_y; j++){
		for(i = 0 ; i < size_x ; i++)
			if(array[size_x * j + i] > 1e-2 || array[size_x * j + i] == 0. ){
				fprintf(fp, " %3.lf", array[size_x * j + i]);
			} else {
				fprintf(fp, " %3.1e", array[size_x * j + i]);
			}
		fprintf(fp, "\n\n");
	}
}

// вывод переменных для дальнейшего построения грайиков в гнуплоте
void print_vars(FILE* f_u13, FILE* f_u14, FILE* f_u11, FILE* f_il,double* vars, double time){
	fprintf(f_u13, "%.10lf %.30lf\n", time, vars[18]);
	fprintf(f_u14, "%.10lf %.30lf\n", time, vars[19]);
	fprintf(f_u11, "%.10lf %.30lf\n", time, vars[20]);
	fprintf(f_il, "%.10lf %.30lf\n", time, vars[6]);
}

void solve(double* mas, double* mas1, double* x){
	x[N - 1] = mas1[N - 1] / mas[(N - 1) * N +  N - 1];
	for (int i = N - 2; i >= 0; i--){
		x[i] = mas1[i];
		for (int j = N - 2; j >= i; j--){
			x[i] = x[i] - mas[i * N + j + 1] * x[j + 1];
		}
		x[i] = x[i] / mas[i * N + i];
	}    
}

void gauss(double* Y, double* B){	   
	double c = 0;

	for (int k = 0; k < N; k++){
		for (int i = k + 1; i < N; i++){
			c = -Y[i * N + k] / Y[k * N + k];
			B[i] = B[i] + B[k] * c;

			for (int j = k; j < N; j++){  			
				Y[i * N + j] = Y[i * N + j] + Y[k * N + j] * c;
				if (k == 19 && i ==21 && j == 19){
					Y[i * N + j] = 0.0;
				}
			}
		}
	}
}

void fill_yakobi_matrix (double* jakobi, double d_dt, double* vars_last_iter){
	int i, j, k, h;
	int diagonal = BRANCHES + HORDES + M;
	double a = - It_d11 * (exp(vars_last_iter[4] / MFt_d11 ) ) / MFt_d11;
	
	init_array(jakobi, N, N);

	for ( k = 0; k < diagonal; k++) {
		for (h = 0; h < diagonal; h++) {
			if(k == h){
				jakobi[N * k + h] = 1.;
			}
		}
	}

	jakobi[N * 0 + 18] = -1./d_dt;
	jakobi[N * 1 + 19] = -1./d_dt;
	jakobi[N * 2 + 20] = -1./d_dt;
	jakobi[N * 3 + 15] = -1./d_dt;

	// матрица M 
	jakobi[N * 4 + 20] = -1.;
	jakobi[N * 5 + 19] = -1.;
	jakobi[N * 5 + 20] = 1.;
	jakobi[N * 5 + 21] = 1.;
	jakobi[N * 6 + 17] = -1.;
	jakobi[N * 6 + 18] = -1.;
	jakobi[N * 6 + 19] = 1.;
	jakobi[N * 7 + 20] = -1.;

	// транспонированная матрица M
	jakobi[N * 8 + 15]  = 1.;
	jakobi[N * 9 + 15]  = 1.;
	jakobi[N * 10 + 14] = 1.;
	jakobi[N * 10 + 15] = -1.;
	jakobi[N * 11 + 13] = 1.;
	jakobi[N * 11 + 14] = -1.;
	jakobi[N * 11 + 16] = 1.;
	jakobi[N * 12 + 14] = -1.;

	// нижняя часть матрицы
	jakobi[N * 13 + 4]  = -1.;
	jakobi[N * 13 + 13] = R_d11;
	jakobi[N * 14 + 5]  = -1.;
	jakobi[N * 14 + 14] = R;
	jakobi[N * 15 + 3]  = -L;
	jakobi[N * 15 + 6]  = 1.;
	jakobi[N * 16 + 4]  = a;
	jakobi[N * 16 + 16] = 1.;
	jakobi[N * 17 + 17] = 1.;
	jakobi[N * 18 + 0]  = -C13;
	jakobi[N * 18 + 9]  = 1.;
	jakobi[N * 19 + 1]  = -C14;
	jakobi[N * 19 + 10] = 1.;
	jakobi[N * 20 + 2]  = -C_d11;
	jakobi[N * 20 + 11] = 1.;
	jakobi[N * 21 + 12] = -Rb_d11;
	jakobi[N * 21 + 21] = 1.;
}

// функция заполняет вектор значений функции
void fill_func_values(double* func_values, double* bVars, double* last_altVars, double dt){
	func_values[0]  = - (bVars[0]  - (bVars[18] - last_altVars[18]) / dt);		 
	func_values[1]  = - (bVars[1]  - (bVars[19] - last_altVars[19]) / dt);
	func_values[2]  = - (bVars[2]  - (bVars[20] - last_altVars[20]) / dt);
	func_values[3]  = - (bVars[3]  - (bVars[15] - last_altVars[15]) / dt);

	func_values[4]	= - (bVars[4]  - bVars[20]);
	func_values[5]  = - (bVars[5]  - bVars[19] + bVars[20] + bVars[21]);
	func_values[6]  = - (bVars[6]  - bVars[17] - bVars[18] + bVars[19]);
	func_values[7]  = - (bVars[7]  - bVars[20]);

	func_values[8]  = - (bVars[8]  + bVars[15]);
	func_values[9]  = - (bVars[9]  + bVars[15]);
	func_values[10] = - (bVars[10] + bVars[14] - bVars[15]);
	func_values[11] = - (bVars[11] + bVars[13] - bVars[14] + bVars[16]);
	func_values[12] = - (bVars[12] - bVars[14]);

	func_values[13] = - (bVars[13] * R_d11 - bVars[4]);
	func_values[14] = - (bVars[14] * R - bVars[5]);
	func_values[15] = - (bVars[6]  - L * bVars[3]);
	func_values[16] = - (bVars[16] - It_d11 * (exp(bVars[4] / MFt_d11) - 1.));
	func_values[17] = - (bVars[17] - E);
	func_values[18] = - (bVars[9]  - C13 * bVars[0] );
	func_values[19] = - (bVars[10] - C14 * bVars[1]);
	func_values[20] = - (bVars[11] - C_d11 * bVars[2]);
	func_values[21] = - (bVars[21] - Rb_d11 * bVars[12]);
}

// заполнение массива базисных переменных 
void fill_base_vars(double* deltas, double* bVars){	
	for(int i = 0; i < N; i++){
		bVars[i] += deltas[i]; 
	}
}

// возвращаемся к изначальным значениям при неудачном шаге
void to_last_step(double* bVars, double* last_step, double* last_but_one){
	for (int i = 0; i < N; i++){
		bVars[i] = last_step[i];
		last_step[i] = last_but_one[i]; 
	}
}

// проверка локальной точности. Если шаг удачный вернет true
int check_local_eps(double* dt, double* last_dt, double* last_step, double* last_but_one, double* current_step){
	double eps1  = 0.001;
	double eps2  = 0.005;      
	double fcur = current_step[4];
	double fprev = last_but_one[4];
	double flast = last_step[4];
	double local_eps = fabs( (*dt / (*dt + *last_dt)) * (fcur - flast - (*dt / *last_dt) * (flast - fprev)));
	double dt_max = 1e-5;

	if (local_eps < eps1){
		(*last_dt) = (*dt);
		(*dt) = 2 * (*dt);
		return true;
	} else if (local_eps < eps2){
		(*last_dt) = (*dt);
		return true;
	} else {
		(*dt) = (*dt) / 2;
	}	
	return false;
}

// сохраняем значения на предыдущем шаге, чтобы можно было проверить точность и вернуться к этим значениям
// при неудачном шаге
void save_last_step(double* current, double* last_step, double* last_but_one){
	for (int i = 0 ; i < N ; i++){
		last_but_one[i] = last_step[i];
		last_step[i] = current[i];
	}	
}

int main(){
	double bVars[N];				// массив, содержащий базисные переменные
	double delta_bVars[N];			// массив, сожержащий дельты базисных переменных
	double func_values[N];			// вектор значений функции
	double jakobi[N * N];			// матрица якоби	
	double last_step_bVars[N];  	// массив, хранящий последний шаг для каждой переменной
	double last_but_one_bVars[N];  	// массив, хранящий предпоследний шаг для каждой переменной
	double dt = 1e-9;				// шаг интегрирования
	double last_dt = 1e-9;			// переменная для хранения предыдущего 
	double max_time = 1e-2;			// максимальное время интегрирования
	double eps = 1e-5;				// eps для метода Ньютона
	double time = 0;
	int endFlag = false;	
	int iteration_num = 0;
	int count_steps = 0;		

	FILE* f_uC13;
	FILE* f_uC14;
	FILE* f_uC11; 
	FILE* f_iL;	
	FILE* f_dt;
	FILE* f_matrix;
	FILE* f_matrix2;		

	f_uC13 = fopen("out/f_uC13.txt", "w");
	f_uC14 = fopen("out/f_uC14.txt", "w");
	f_uC11 = fopen("out/f_uC11.txt", "w");
	f_iL   = fopen("out/f_iL.txt", "w");
	f_dt   = fopen("out/f_dt.txt", "w");
	f_matrix = fopen("out/f_matrix.txt", "w");
	f_matrix2 = fopen("out/f_matrix2.txt", "w");
	
	init_array(bVars, 1, N);
	init_array(delta_bVars, 1, N);	
	init_array(last_step_bVars, 1, N);	
	init_array(last_but_one_bVars, 1, N);

	while(time < max_time) {
		
		iteration_num = 0;

		do {
			fill_yakobi_matrix(jakobi, dt, bVars);
			fill_func_values(func_values, bVars, last_step_bVars, dt);						
			gauss(jakobi, func_values);

			if(count_steps < 100) {				
				fprintf(f_matrix,"time = %.12e, iter = %d, dt = %e\n", time, iteration_num, dt);
				print_array_to_file(jakobi, N, N,f_matrix);
				fprintf(f_matrix, "Deltas: \t Functions\n");
				for(int j = 0; j < N; j++){
					fprintf(f_matrix, "%e\t%e\n", delta_bVars[j], func_values[j]);
				}
			} 

			solve(jakobi, func_values, delta_bVars);										
			fill_base_vars(delta_bVars, bVars);
          
			iteration_num += 1;

			endFlag = true;
			for(int i = 4; i < N; i++){
				if(fabs(delta_bVars[i]) >= eps){
					endFlag = false;
				}
			}

			if (endFlag){
				endFlag = check_local_eps(&dt, &last_dt, last_step_bVars, last_but_one_bVars, bVars);				
			}

			if(!endFlag){
				if(iteration_num >= 6){
					dt /= 2;
					to_last_step(bVars, last_step_bVars, last_but_one_bVars);
					iteration_num = 0;
				}
			}
		} while(!endFlag);

		// Соответсвенно тут уже логика сохранения успешного шага.
		/*endFlag = check_local_eps(&dt, &last_dt, last_step_bVars, last_but_one_bVars, bVars);
		if(endFlag){
			// шаг удался 
			// отладочный вывод
			if(count_steps % 10000 == 0){
				printf("Time:%e, Step:%d,  dt = %e\n", time ,count_steps, dt );												
			}
 
			save_last_step(bVars, last_step_bVars,  last_but_one_bVars);			
			print_vars(f_uC13, f_uC14, f_uC11, f_iL, bVars, time);
			time += dt;
			count_steps += 1;			
		}
		else{			
			// шаг был неудачным
			to_last_step(bVars, last_step_bVars,  last_but_one_bVars);
		}
		*/
		if(count_steps % 10000 == 0){
			printf("Time:%e, Step:%d,  dt = %e\n", time ,count_steps, dt );
		}

		save_last_step(bVars, last_step_bVars, last_but_one_bVars);
		print_vars(f_uC13, f_uC14, f_uC11, f_iL, bVars, time);
		time += dt;
		count_steps += 1;
	}

	fclose(f_uC13);
	fclose(f_uC14);
	fclose(f_uC11); 
	fclose(f_iL);	
	fclose(f_dt);
	fclose(f_matrix);	
	fclose(f_matrix2);	
	return 0;
}