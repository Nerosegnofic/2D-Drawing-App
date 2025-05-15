#include <windows.h>
#include <vector>
#include <algorithm>
#include <array>

using std::array, std::clamp, std::vector;

class Point {
public:
    float x, y;
    COLORREF color;

    Point() : x(0), y(0), color(RGB(0, 0, 0)) {}
    Point(const float _x, const float _y) : x(_x), y(_y), color(RGB(0, 0, 0)) {}
    Point(const float _x, const float _y, const COLORREF _color) : x(_x), y(_y), color(_color) {}

    Point operator*(const float scalar) const {
        return Point(x * scalar, y * scalar, color);
    }

    Point operator/(const float scalar) const {
        if (scalar == 0.0f) {
            return *this;
        }
        return *this * (1.0f / scalar);
    }

    Point operator+(const Point& p) const {
        return Point(x + p.x, y + p.y);
    }

    Point operator-(const Point& p) const {
        return Point(x - p.x, y - p.y);
    }
};

class BezierCurve {
    static constexpr array<array<int, 4>, 4> BEZIER_MATRIX {{
        {-1,  3, -3, 1},
        { 3, -6,  3, 0},
        {-3,  3,  0, 0},
        { 1,  0,  0, 0}
    }};

    static int Round(const float value) {
        return static_cast<int>(value + 0.5f);
    }

    static vector<vector<int>> matrix_multiply(
        const array<array<int, 4>, 4>& A,
        const vector<vector<int>>& B
    ) {
        const size_t rows {A.size()};
        const size_t cols {B[0].size()};
        const size_t inner {B.size()};

        vector result(rows, vector(cols, 0));

        for (size_t i {0}; i < rows; ++i) {
            for (size_t j {0}; j < cols; ++j) {
                for (size_t k {0}; k < inner; ++k) {
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return result;
    }

public:
    static void draw(const HDC hdc, const std::vector<Point>& control_points, const int number_of_segments) {
        if (control_points.size() != 4) {
            return;
        }

        vector values(4, vector<int>(5));
        for (int i {0}; i < 4; ++i) {
            values[i][0] = Round(control_points[i].x);
            values[i][1] = Round(control_points[i].y);
            values[i][2] = GetRValue(control_points[i].color);
            values[i][3] = GetGValue(control_points[i].color);
            values[i][4] = GetBValue(control_points[i].color);
        }

        const auto coefficients {matrix_multiply(BEZIER_MATRIX, values)};

        const float step {1.0f / number_of_segments};
        for (float t {0.0f}; t <= 1.0f; t += step) {
            const float t2 {t * t};
            const float t3 {t2 * t};

            vector result(5, 0.0f);
            for (int i {0}; i < 5; ++i) {
                result[i] = t3 * coefficients[0][i] +
                            t2 * coefficients[1][i] +
                            t  * coefficients[2][i] +
                                 coefficients[3][i];
            }

            SetPixel(
                hdc,
                Round(result[0]),
                Round(result[1]),
                RGB(
                    clamp(Round(result[2]), 0, 255),
                    clamp(Round(result[3]), 0, 255),
                    clamp(Round(result[4]), 0, 255)
                )
            );
        }
    }
};

class BezierApp {
    static constexpr int MAX_CONTROL_POINTS {4};
    static constexpr int DEFAULT_CURVE_SEGMENTS {2000};

    static inline vector<Point> control_points;
    static inline int point_count {0};

    static constexpr COLORREF CONTROL_POINT_COLORS[MAX_CONTROL_POINTS] = {
        RGB(0, 0, 0),
        RGB(255, 0, 0),
        RGB(0, 0, 255),
        RGB(0, 255, 0)
    };

public:
    static void on_click(const HWND hwnd, const int x, const int y) {
        if (point_count < MAX_CONTROL_POINTS) {
            control_points.push_back(Point(
                static_cast<float>(x),
                static_cast<float>(y),
                CONTROL_POINT_COLORS[point_count]
            ));

            ++point_count;

            if (point_count == MAX_CONTROL_POINTS) {
                const HDC hdc {GetDC(hwnd)};
                BezierCurve::draw(hdc, control_points, DEFAULT_CURVE_SEGMENTS);
                ReleaseDC(hwnd, hdc);

                control_points.clear();
                point_count = 0;
            }
        }
    }
};

LRESULT CALLBACK WndProc(const HWND hwnd, const UINT msg, const WPARAM wParam, const LPARAM lParam) {
    switch (msg) {
        case WM_LBUTTONDOWN:
            BezierApp::on_click(
                hwnd,
                LOWORD(lParam),
                HIWORD(lParam)
            );
            break;

        case WM_DESTROY:
            PostQuitMessage(0);
            break;

        default:
            return DefWindowProc(hwnd, msg, wParam, lParam);
    }
    return 0;
}

int WINAPI WinMain(const HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, const int nShowCmd) {
    WNDCLASSW wc {};
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = L"BezierCurveApp";
    wc.hbrBackground = reinterpret_cast<HBRUSH>(COLOR_WINDOW + 1);
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);

    if (!RegisterClassW(&wc)) {
        MessageBoxW(nullptr, L"Window Registration Failed!", L"Error", MB_ICONEXCLAMATION | MB_OK);
        return 0;
    }

    const HWND hwnd {CreateWindowW(
        L"BezierCurveApp",
        L"Bezier Curve Drawing App",
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, 800, 600,
        NULL, NULL, hInstance, NULL
    )};

    if (hwnd == nullptr) {
        MessageBoxW(nullptr, L"Window Creation Failed!", L"Error", MB_ICONEXCLAMATION | MB_OK);
        return 0;
    }

    ShowWindow(hwnd, nShowCmd);
    UpdateWindow(hwnd);

    MSG msg {};
    while (GetMessage(&msg, nullptr, 0, 0)) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    return static_cast<int>(msg.wParam);
}
