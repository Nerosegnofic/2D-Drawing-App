#include <math.h>
#include <Windows.h>

typedef struct {
    int x, y;
} Point;

Point points[5];
int point_count = 0;

void draw_line(const HDC hdc, const int x1, const int y1, const int x2, const int y2, const COLORREF c) {
    const int dx = x2 - x1;
    const int dy = y2 - y1;

    const int steps = abs(dx) > abs(dy) ? abs(dx) : abs(dy);

    const float x_inc = dx / (float)steps;
    const float y_inc = dy / (float)steps;

    float x = x1;
    float y = y1;

    for (int i = 0; i <= steps; ++i) {
        SetPixel(hdc, round(x), round(y), c);
        x += x_inc;
        y += y_inc;
    }
}

LRESULT CALLBACK WndProc(const HWND hwnd, const UINT m, const WPARAM wp, const LPARAM lp) {
    switch (m) {
        case WM_LBUTTONDOWN:
            if (point_count < 5) {
                points[point_count].x = LOWORD(lp);
                points[point_count].y = HIWORD(lp);
                ++point_count;
            }
            if (point_count == 5) {
                const HDC hdc = GetDC(hwnd);
                draw_line(hdc, points[0].x, points[0].y, points[2].x, points[2].y, RGB(255, 0, 0));
                draw_line(hdc, points[1].x, points[1].y, points[3].x, points[3].y, RGB(0, 255, 0));
                draw_line(hdc, points[2].x, points[2].y, points[4].x, points[4].y, RGB(0, 0, 255));
                draw_line(hdc, points[3].x, points[3].y, points[0].x, points[0].y, RGB(0, 0, 0));
                draw_line(hdc, points[4].x, points[4].y, points[1].x, points[1].y, RGB(139, 0, 139));
                ReleaseDC(hwnd, hdc);
                point_count = 0;
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
    wc.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    wc.hIcon = LoadIcon(nullptr, IDI_APPLICATION);
    wc.lpszClassName = L"MyClass";
    wc.lpfnWndProc = WndProc;
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.hInstance = hInstance;
    RegisterClass(&wc);

    const HWND hwnd = CreateWindow(L"MyClass", L"Draw 5-Point Lines", WS_OVERLAPPEDWINDOW,
                             CW_USEDEFAULT, CW_USEDEFAULT, 500, 500,
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
