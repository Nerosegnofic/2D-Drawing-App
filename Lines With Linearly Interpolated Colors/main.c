#include <Windows.h>

int round(const double x) {
    return (int)(x + 0.5);
}

void interpolated_colored_line(const HDC hdc, const int x1, const int y1, const int x2, const int y2, const COLORREF c1, const COLORREF c2) {
    const int dx = x2 - x1;
    const int dy = y2 - y1;
    const int steps = max(abs(dx), abs(dy));

    const double r1 = GetRValue(c1);
    const double g1 = GetGValue(c1);
    const double b1 = GetBValue(c1);
    const double r2 = GetRValue(c2);
    const double g2 = GetGValue(c2);
    const double b2 = GetBValue(c2);

    const double dr = (r2 - r1) / steps;
    const double dg = (g2 - g1) / steps;
    const double db = (b2 - b1) / steps;

    double x = x1, y = y1;
    double r = r1, g = g1, b = b1;

    for (int i = 0; i <= steps; ++i) {
        SetPixel(hdc, round(x), round(y), RGB(round(r), round(g), round(b)));
        x += (double)dx / steps;
        y += (double)dy / steps;
        r += dr;
        g += dg;
        b += db;
    }
}

POINT start_point;
bool is_drawing = false;

LRESULT CALLBACK WndProc(const HWND hwnd, const UINT m, const WPARAM wp, const LPARAM lp) {
    switch (m) {
        case WM_LBUTTONDOWN:
            start_point.x = LOWORD(lp);
            start_point.y = HIWORD(lp);
            is_drawing = true;
            break;

        case WM_LBUTTONUP:
            if (is_drawing) {
                const HDC hdc = GetDC(hwnd);
                const int x2 = LOWORD(lp);
                const int y2 = HIWORD(lp);
                interpolated_colored_line(hdc, start_point.x, start_point.y, x2, y2, RGB(255, 0, 0), RGB(0, 0, 255));
                ReleaseDC(hwnd, hdc);
                is_drawing = false;
            }
            break;

        case WM_CLOSE:
            DestroyWindow(hwnd);
            break;

        case WM_DESTROY:
            PostQuitMessage(0);
            break;

        default:
            return DefWindowProc(hwnd, m, wp, lp);
    }
    return 0;
}

int APIENTRY WinMain(const HINSTANCE hInstance, [[maybe_unused]] HINSTANCE hPrevInstance, [[maybe_unused]] LPSTR lpCmdLine, const int nShowCmd) {
    WNDCLASS wc = {};
    wc.hbrBackground = (HBRUSH)GetStockObject(LTGRAY_BRUSH);
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.hIcon = LoadIcon(nullptr, IDI_WINLOGO);
    wc.lpszClassName = L"MyClass";
    wc.lpfnWndProc = WndProc;
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.hInstance = hInstance;
    RegisterClass(&wc);

    const HWND hwnd = CreateWindow(L"MyClass", L"Draw Line", WS_OVERLAPPEDWINDOW,
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
