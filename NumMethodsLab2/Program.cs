using System;

class NumericalMethods
{
    static void Main()
    {
        double[,] matrixGauss = {
            { 4, 0, 1, 0 },
            { 0, 3, 0, 2 },
            { 1, 0, 5, 1 },
            { 0, 2, 1, 4 }
        };
        double[] vectorGauss = { 12, 19, 27, 30 };

        Console.WriteLine("The system of equations for the Gaussian method:");
        PrintSystem(matrixGauss, vectorGauss);

        SolveGauss(matrixGauss, vectorGauss);

        double[,] matrixCholesky = {
            { 1, 2, 0 },
            { 2, 2, 4 },
            { 0, 4, 3 }
        };
        double[] vectorCholesky = { 5, 22, 20 };

        Console.WriteLine("\nThe system of equations for the Cholesky method:");
        PrintSystem(matrixCholesky, vectorCholesky);

        SolveCholesky(matrixCholesky, vectorCholesky);

        double[,] matrixSeidel = {
            { 5, 1, 1, 0 },
            { 1, 2, 0, 0 },
            { 1, 0, 4, 2 },
            { 0, 0, 2, 3 }
        };
        double[] vectorSeidel = { 17, 8, 28, 23 };

        Console.WriteLine("\nThe system of equations for the Seidel method:");
        PrintSystem(matrixSeidel, vectorSeidel);

        SolveSeidel(matrixSeidel, vectorSeidel);
    }

    static void PrintSystem(double[,] matrix, double[] vector)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                Console.Write($"{matrix[i, j]} ");
            }
            Console.WriteLine($"= {vector[i]}");
        }
    }

    static void SolveGauss(double[,] matrix, double[] vector)
    {
        int n = vector.Length;
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                double factor = matrix[j, i] / matrix[i, i];
                for (int k = i; k < n; k++)
                    matrix[j, k] -= factor * matrix[i, k];
                vector[j] -= factor * vector[i];
            }
        }

        double[] result = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            result[i] = vector[i];
            for (int j = i + 1; j < n; j++)
                result[i] -= matrix[i, j] * result[j];
            result[i] /= matrix[i, i];
        }

        Console.WriteLine("\nSolution by the Gaussian method:");
        for (int i = 0; i < n; i++)
        {
            Console.WriteLine($"X{i + 1} = {result[i]:F4}");
        }
    }

    static void SolveCholesky(double[,] matrix, double[] vector)
    {
        int n = vector.Length;
        double[,] L = new double[n, n];
        double[] result = new double[n];

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                double sum = 0;
                for (int k = 0; k < j; k++)
                    sum += L[i, k] * L[j, k];
                if (i == j)
                    L[i, j] = Math.Sqrt(matrix[i, i] - sum);
                else
                    L[i, j] = (matrix[i, j] - sum) / L[j, j];
            }
        }

        for (int i = 0; i < n; i++)
        {
            result[i] = vector[i];
            for (int j = 0; j < i; j++)
                result[i] -= L[i, j] * result[j];
            result[i] /= L[i, i];
        }

        for (int i = n - 1; i >= 0; i--)
        {
            for (int j = i + 1; j < n; j++)
                result[i] -= L[j, i] * result[j];
            result[i] /= L[i, i];
        }

        Console.WriteLine("\nSolution by the Cholesky method:");
        for (int i = 0; i < n; i++)
        {
            Console.WriteLine($"X{i + 1} = {result[i]:F4}");
        }
    }

    static void SolveSeidel(double[,] matrix, double[] vector)
    {
        int n = vector.Length;
        double[] x = new double[n];
        double[] prevX = new double[n];
        double tolerance = 1e-6; 
        bool stop = false;

        while (!stop)
        {
            stop = true;
            for (int i = 0; i < n; i++)
            {
                double sum = vector[i];
                for (int j = 0; j < n; j++)
                {
                    if (i != j)
                        sum -= matrix[i, j] * x[j];
                }
                prevX[i] = x[i];
                x[i] = sum / matrix[i, i];

                if (Math.Abs(x[i] - prevX[i]) > tolerance)
                    stop = false;
            }
        }

        Console.WriteLine("\nSolution by the Seidel method:");
        for (int i = 0; i < n; i++)
        {
            Console.WriteLine($"X{i + 1} = {x[i]:F4}");
        }
    }
}
