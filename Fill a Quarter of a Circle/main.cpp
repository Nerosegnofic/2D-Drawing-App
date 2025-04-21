#define UNICODE
#define UNICODE
#include <Windows.h>
#include <cmath>

bool is_first_circle {true};

void draw_line(const HDC hdc, const int x1, const int y1, const int x2, const int y2, const COLORREF c) {
    const int dx {x2 - x1};
    const int dy {y2 - y1};

    const int steps {abs(dx) > abs(dy) ? abs(dx) : abs(dy)};

    const float x_inc {dx / static_cast<float>(steps)};
    const float y_inc {dy / static_cast<float>(steps)};

    float x = x1;
    float y = y1;

    for (int i {0}; i <= steps; ++i) {
        SetPixel(hdc, std::round(x), std::round(y), c);
        x += x_inc;
        y += y_inc;
    }
}

void draw_circle_points(const HDC hdc, const int xc, const int yc, const int a, const int b, const COLORREF c) {
    SetPixel(hdc, xc + a, yc + b, c);
    SetPixel(hdc, xc - a, yc + b, c);
    SetPixel(hdc, xc + a, yc - b, c);
    SetPixel(hdc, xc - a, yc - b, c);

    SetPixel(hdc, xc + b, yc + a, c);
    SetPixel(hdc, xc - b, yc + a, c);
    SetPixel(hdc, xc + b, yc - a, c);
    SetPixel(hdc, xc - b, yc - a, c);

    if (is_first_circle) {
        if (a >= 0 && b >= 0) {
            draw_line(hdc, xc, yc, xc + a, yc - b, RGB(0, 0, 0));
            draw_line(hdc, xc, yc, xc + b, yc - a, RGB(0, 0, 0));
        }
    }
}

void draw_bresenham_circle(const HDC hdc, const int xc, const int yc, const int r, const COLORREF c) {
    int x {0};
    int y {r};
    int d {1 - r};
    int d1 {3};
    int d2 {5 - 2 * r};
    draw_circle_points(hdc, xc, yc, x, y, c);

    while (x < y) {
        if (d < 0) {
            d += d1;
            d2 += 2;
        } else {
            d += d2;
            d2 += 4;
            --y;
        }
        d1 += 2;
        ++x;
        draw_circle_points(hdc, xc, yc, x, y, c);
    }
}

LRESULT WndProc(const HWND hwnd, const UINT msg, const WPARAM wp, const LPARAM lp) {
    static POINT center1 {};
    static POINT center2 {};
    static int click_count {0};

    switch (msg) {
    case WM_LBUTTONDOWN:
        if (click_count == 0) {
            center1.x = LOWORD(lp);
            center1.y = HIWORD(lp);
            ++click_count;
        } else if (click_count == 1) {
            center2.x = LOWORD(lp);
            center2.y = HIWORD(lp);
            ++click_count;
        }
        break;

    case WM_RBUTTONDOWN:
        if (click_count == 2) {
            POINT perimeter;
            perimeter.x = LOWORD(lp);
            perimeter.y = HIWORD(lp);

            const int dx1 {static_cast<int>(center1.x - perimeter.x)};
            const int dy1 {static_cast<int>(center1.y - perimeter.y)};
            const int r1 {static_cast<int>(round(sqrt(dx1 * dx1 + dy1 * dy1)))};

            const int dx2 {static_cast<int>(center2.x - perimeter.x)};
            const int dy2 {static_cast<int>(center2.y - perimeter.y)};
            const int r2 {static_cast<int>(round(sqrt(dx2 * dx2 + dy2 * dy2)))};

            const HDC hdc {GetDC(hwnd)};
            draw_bresenham_circle(hdc, center1.x, center1.y, r1, RGB(0, 0, 0));
            is_first_circle = false;
            draw_bresenham_circle(hdc, center2.x, center2.y, r2, RGB(0, 0, 0));
            ReleaseDC(hwnd, hdc);

            click_count = 0;
            is_first_circle = true;
        }
        break;

    case WM_CLOSE:
        DestroyWindow(hwnd);
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hwnd, msg, wp, lp);
    }
    return 0;
}

int APIENTRY WinMain(const HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, const int nShowCmd) {
    WNDCLASS wc {};
    wc.cbClsExtra = 0;
    wc.cbWndExtra = 0;
    wc.hbrBackground = static_cast<HBRUSH>(GetStockObject(LTGRAY_BRUSH));
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.hIcon = LoadIcon(nullptr, IDI_WINLOGO);
    wc.lpszClassName = L"MyClass";
    wc.lpszMenuName = nullptr;
    wc.lpfnWndProc = WndProc;
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.hInstance = hInstance;
    RegisterClass(&wc);

    HWND hwnd = CreateWindow(L"MyClass", L"Dual Circle Drawer", WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT,
        nullptr, nullptr, hInstance, nullptr);
    ShowWindow(hwnd, nShowCmd);
    UpdateWindow(hwnd);

    MSG msg;
    while (GetMessage(&msg, nullptr, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    return 0;
}
