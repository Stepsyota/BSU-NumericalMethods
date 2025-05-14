#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

const double EPSILON1 = 1e-3;
const double EPSILON2 = 1e-5;
const int MAX_ITERATIONS = 1000;


void printHorizontalLine() {
    cout << "========================================================" << endl;
}

// Явный метод Эйлера для решения системы ОДУ
void solveExplicitEulerSystem(double parameter_a, double parameter_omega, double epsilon = EPSILON1) {
    const double START_TIME = 0.0;
    const double END_TIME = 1.0;
    const double MAX_STEP = 0.1;

    vector<double> solution = { 0.0, -0.412 }; // Начальные условия
    double current_time = START_TIME;

    printHorizontalLine();
    cout << "Решение системы ОДУ явным методом Эйлера\n";
    cout << "Параметр a = " << parameter_a << ", частота w = " << parameter_omega << "\n";
    printHorizontalLine();

    cout << setw(10) << "Время" << setw(15) << "u1" << setw(15) << "u2" << setw(15) << "Шаг" << endl;
    cout << fixed << setprecision(6);

    while (current_time < END_TIME) {
        // Вычисление правых частей уравнений
        double right_part1 = -solution[0] * solution[1] + (current_time == 0 ? 1.0 : sin(current_time) / current_time);
        double right_part2 = -solution[1] * solution[1] + (parameter_a * current_time) / (1 + current_time * current_time);

        // Вычисление допустимых шагов
        double max_u1 = max(1.0, abs(solution[0]));
        double max_u2 = max(1.0, abs(solution[1]));
        double local_epsilon1 = 0.01 * max_u1;
        double local_epsilon2 = 0.01 * max_u2;

        double step1 = (right_part1 != 0) ? local_epsilon1 / abs(right_part1) : MAX_STEP;
        double step2 = (right_part2 != 0) ? local_epsilon2 / abs(right_part2) : MAX_STEP;
        double current_step = min({ step1, step2, MAX_STEP });

        // Коррекция шага при выходе за конечное время
        if (current_time + current_step > END_TIME) current_step = END_TIME - current_time;

        // Выполнение шага метода
        solution[0] += current_step * right_part1;
        solution[1] += current_step * right_part2;
        current_time += current_step;

        // Вывод результатов
        cout << setw(10) << current_time << setw(15) << solution[0]
            << setw(15) << solution[1] << setw(15) << current_step << endl;
    }
}

// Функция для вычисления правых частей системы (матрица A и вектор b)
vector<double> computeSystemRightParts(const vector<double>& solution, const vector<double>& eigenvalues) {
    double lambda1 = eigenvalues[0], lambda2 = eigenvalues[1], lambda3 = eigenvalues[2];

    // Матрица коэффициентов A
    vector<vector<double>> matrix_A = {
        {(2 * lambda1 + 4 * lambda2) / 6, (2 * (lambda1 - lambda2)) / 6, (2 * (lambda1 - lambda2)) / 6},
        {(2 * (lambda1 - lambda2)) / 6, (2 * lambda1 + lambda2 + 3 * lambda3) / 6, (2 * lambda1 + lambda2 - 3 * lambda3) / 6},
        {(2 * (lambda1 - lambda2)) / 6, (2 * lambda1 + lambda2 - 3 * lambda3) / 6, (2 * lambda1 + lambda2 + 3 * lambda3) / 6}
    };

    // Вектор правой части b
    vector<double> vector_b = {
        -(4 * lambda1 + 2 * lambda2) / 6,
        -(4 * lambda1 - lambda2 - 9 * lambda3) / 6,
        -(4 * lambda1 - lambda2 + 9 * lambda3) / 6
    };

    // Вычисление f(u) = A*u - b
    vector<double> result(3);
    for (int i = 0; i < 3; ++i) {
        result[i] = -vector_b[i];
        for (int j = 0; j < 3; ++j) {
            result[i] += matrix_A[i][j] * solution[j];
        }
    }
    return result;
}

// Метод Ньютона для решения нелинейной системы
vector<double> solveByNewtonMethod(const vector<double>& prev_solution, const vector<double>& eigenvalues,
    double step_size, double epsilon = 1e-6, int max_iterations = 100) {
    vector<double> current_solution = prev_solution;
    int system_size = current_solution.size();

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        vector<double> F = computeSystemRightParts(current_solution, eigenvalues);

        // Вычисление невязки F(u) = u - u_prev - step_size*f(u)
        for (int i = 0; i < system_size; ++i) {
            F[i] = current_solution[i] - prev_solution[i] - step_size * F[i];
        }

        // Проверка условия сходимости
        double residual_norm = 0.0;
        for (double value : F) residual_norm += value * value;
        residual_norm = sqrt(residual_norm);
        if (residual_norm < epsilon) break;

        // Построение матрицы Якоби J = I - step_size*df/du
        vector<vector<double>> Jacobian(system_size, vector<double>(system_size, 0.0));
        for (int i = 0; i < system_size; ++i) Jacobian[i][i] = 1.0;

        double lambda1 = eigenvalues[0], lambda2 = eigenvalues[1], lambda3 = eigenvalues[2];
        vector<vector<double>> matrix_A = {
            {(2 * lambda1 + 4 * lambda2) / 6, (2 * (lambda1 - lambda2)) / 6, (2 * (lambda1 - lambda2)) / 6},
            {(2 * (lambda1 - lambda2)) / 6, (2 * lambda1 + lambda2 + 3 * lambda3) / 6, (2 * lambda1 + lambda2 - 3 * lambda3) / 6},
            {(2 * (lambda1 - lambda2)) / 6, (2 * lambda1 + lambda2 - 3 * lambda3) / 6, (2 * lambda1 + lambda2 + 3 * lambda3) / 6}
        };

        for (int i = 0; i < system_size; ++i) {
            for (int j = 0; j < system_size; ++j) {
                Jacobian[i][j] -= step_size * matrix_A[i][j];
            }
        }

        // Решение системы J*delta = -F методом Гаусса
        vector<double> delta = F;
        for (int i = 0; i < system_size; ++i) delta[i] = -delta[i];

        // Прямой ход метода Гаусса
        for (int i = 0; i < system_size; ++i) {
            // Поиск ведущего элемента
            int pivot_row = i;
            for (int k = i + 1; k < system_size; ++k) {
                if (abs(Jacobian[k][i]) > abs(Jacobian[pivot_row][i])) {
                    pivot_row = k;
                }
            }

            // Перестановка строк
            swap(Jacobian[i], Jacobian[pivot_row]);
            swap(delta[i], delta[pivot_row]);

            // Исключение переменных
            for (int k = i + 1; k < system_size; ++k) {
                double factor = Jacobian[k][i] / Jacobian[i][i];
                for (int j = i; j < system_size; ++j) {
                    Jacobian[k][j] -= factor * Jacobian[i][j];
                }
                delta[k] -= factor * delta[i];
            }
        }

        // Обратный ход метода Гаусса
        for (int i = system_size - 1; i >= 0; --i) {
            for (int j = i + 1; j < system_size; ++j) {
                delta[i] -= Jacobian[i][j] * delta[j];
            }
            delta[i] /= Jacobian[i][i];
        }

        // Обновление решения
        for (int i = 0; i < system_size; ++i) {
            current_solution[i] += delta[i];
        }
    }

    return current_solution;
}

// Неявный метод Эйлера с выбором шага по стратегии
void solveImplicitEulerSystem(const vector<double>& eigenvalues, double end_time,
    double epsilon = EPSILON1, bool use_three_zone_strategy = true) {
    const double START_TIME = 0.0;
    const double MIN_STEP = 1e-6; // Минимальный шаг уменьшен для eps2
    const double MAX_STEP = 0.1;
    const int MAX_STEP_REDUCTIONS = 20; // Максимальное количество уменьшений шага подряд

    vector<double> solution = { 10.0, 22.0, 9.0 };
    vector<double> prev_solution = solution;
    vector<double> prev_prev_solution = solution;

    double current_time = START_TIME;
    double prev_step = MIN_STEP;
    double current_step = MIN_STEP;
    int step_reduction_count = 0;
    int iteration_count = 0;

    printHorizontalLine();
    cout << "Решение системы ОДУ неявным методом Эйлера (eps = " << epsilon << ")\n";
    cout << "Собственные значения lambda = (" << eigenvalues[0] << ", "
        << eigenvalues[1] << ", " << eigenvalues[2] << ")\n";
    cout << "Конечное время T = " << end_time << "\n";
    cout << "Стратегия: " << (use_three_zone_strategy ? "три зоны" : "квазиоптимальная") << "\n";
    printHorizontalLine();

    cout << setw(10) << "Время" << setw(15) << "u1" << setw(15) << "u2"
        << setw(15) << "u3" << setw(15) << "Шаг" << endl;
    cout << fixed << setprecision(6);

    while (current_time < end_time && iteration_count < MAX_ITERATIONS) {
        iteration_count++;

        // Коррекция шага при выходе за конечное время
        if (current_time + current_step > end_time) {
            current_step = end_time - current_time;
        }

        // Решение нелинейной системы методом Ньютона
        vector<double> new_solution = solveByNewtonMethod(solution, eigenvalues, current_step);

        // Вычисление оценки погрешности
        vector<double> error_estimate(3);
        for (int i = 0; i < 3; ++i) {
            error_estimate[i] = -current_step / (current_step + prev_step) *
                (new_solution[i] - solution[i] - (current_step / prev_step) *
                    (solution[i] - prev_solution[i]));
        }

        // Вычисление допустимых погрешностей
        vector<double> allowed_error(3);
        for (int i = 0; i < 3; ++i) {
            allowed_error[i] = epsilon * max(1.0, abs(new_solution[i]));
        }

        // Проверка условия точности
        bool need_recalculation = false;
        for (int i = 0; i < 3; ++i) {
            if (abs(error_estimate[i]) > allowed_error[i]) {
                need_recalculation = true;
                break;
            }
        }

        if (need_recalculation) {
            // Уменьшение шага и повтор расчета
            current_step /= 2.0;
            step_reduction_count++;

            if (current_step < MIN_STEP) {
                current_step = MIN_STEP;
                // Если шаг уже минимальный, но точность не достигнута,
                // продолжаем с этим шагом, чтобы избежать бесконечного цикла
                if (step_reduction_count > MAX_STEP_REDUCTIONS) {
                    cout << "Предупреждение: Достигнут минимальный шаг " << MIN_STEP
                        << " без достижения требуемой точности." << endl;
                    need_recalculation = false;
                }
            }

            if (need_recalculation) {
                continue;
            }
        }

        // Сброс счетчика уменьшений шага после успешного шага
        step_reduction_count = 0;

        // Выбор нового шага по выбранной стратегии
        vector<double> new_step_candidates(3);
        for (int i = 0; i < 3; ++i) {
            if (use_three_zone_strategy) {
                // Стратегия трех зон
                if (abs(error_estimate[i]) > allowed_error[i]) {
                    new_step_candidates[i] = current_step / 2.0;
                }
                else if (abs(error_estimate[i]) > allowed_error[i] / 4.0) {
                    new_step_candidates[i] = current_step;
                }
                else {
                    new_step_candidates[i] = 2.0 * current_step;
                }
            }
            else {
                // Квазиоптимальная стратегия
                if (abs(error_estimate[i]) < numeric_limits<double>::min()) {
                    new_step_candidates[i] = MAX_STEP;
                }
                else {
                    new_step_candidates[i] = sqrt(allowed_error[i] / abs(error_estimate[i])) * current_step;
                }
            }
        }

        // Выбор минимального шага из кандидатов
        double next_step = *min_element(new_step_candidates.begin(), new_step_candidates.end());
        next_step = max(MIN_STEP, min(MAX_STEP, next_step));

        // Вывод результатов
        cout << setw(10) << current_time + current_step << setw(15) << new_solution[0]
            << setw(15) << new_solution[1] << setw(15) << new_solution[2]
            << setw(15) << current_step << endl;

        // Обновление переменных для следующего шага
        prev_prev_solution = prev_solution;
        prev_solution = solution;
        solution = new_solution;
        current_time += current_step;
        prev_step = current_step;
        current_step = next_step;
    }

    if (iteration_count >= MAX_ITERATIONS) {
        cout << "Предупреждение: Достигнуто максимальное количество итераций ("
            << MAX_ITERATIONS << "). Расчет остановлен." << endl;
    }
}

int main() {
    setlocale(LC_ALL, "Russian");
    // Решение первой задачи явным методом Эйлера
    double omega = 30;
    double parameter_a = 2.5 + omega / 40.0;
    //solveExplicitEulerSystem(parameter_a, omega, EPSILON1);

    // Решение четвертой задачи неявным методом Эйлера
    vector<double> eigenvalues = { -1.0, -2.0, -3.0 };

    // Сначала с стратегией трех зон
    solveImplicitEulerSystem(eigenvalues, 1.0, EPSILON1, true);
    solveImplicitEulerSystem(eigenvalues, 1.0, EPSILON2, true);
    // Затем с квазиоптимальной стратегией
    solveImplicitEulerSystem(eigenvalues, 1.0, EPSILON1, false);
    solveImplicitEulerSystem(eigenvalues, 1.0, EPSILON2, false);


    return 0;
}