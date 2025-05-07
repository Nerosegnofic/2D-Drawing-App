#include <algorithm>
#include <Windows.h>
#include <cmath>
#include <bits/ranges_algo.h>

int Round(double x)
{
    return (int)(x + 0.5);
}

struct Point {
    float x, y;
    Point(float x = 0, float y = 0) : x(x), y(y) {}
    Point operator*(float f) { return Point(x * f, y * f); }
    Point operator/(int f) { return Point(x / f, y / f); }
    Point operator+(Point p) { return Point(x + p.x, y + p.y); }
    Point operator-(Point p) { return Point(x - p.x, y - p.y); }
};

Point *mul(int H[4][4], Point p[4]) {
    Point *res = new Point[4];
    for (int i = 0; i < 4; i++) {
        res[i] = Point(0, 0);
        for (int j = 0; j < 4; j++) {
            res[i] = res[i] + p[j] * H[i][j];
        }
    }
    return res;
}

float distance(Point p1, Point p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow( p1.y - p2.y, 2));
}

Point recBezier1(float t, Point points[], int si, int ei) {
    if (si == ei) return points[si];
    return recBezier1(t, points, si, ei - 1) * (1 - t) + recBezier1(t, points, si + 1, ei) * t;
}

void DrawBezier1(HDC hdc, Point points[], int n, int numOfSteps, COLORREF c) {
    float step = 1.0 / numOfSteps;
    for (float i = 0; i <= 1; i += step) {
        Point p = recBezier1(i, points, 0, n - 1);
        SetPixel(hdc, Round(p.x), Round(p.y), c);
    }
}

void DrawBezier2(HDC hdc, Point p[], COLORREF c) {
    if (distance(p[0], p[3]) <= 1.0)
        return;

    Point q1 = (p[0] + p[1]) / 2, q2 = (p[1] + p[2]) / 2, q3 = (p[2] + p[3]) / 2;
    Point r1 = (q1 + q2) / 2, r2 = (q2 + q3) / 2;
    Point mid = (r1 + r2) / 2;

    SetPixel(hdc, Round(mid.x), Round(mid.y), c);

    Point p1[] = {p[0], q1, r1, mid}, p2[] = {mid, r2, q3, p[3]};

    DrawBezier2(hdc, p1, c);
    DrawBezier2(hdc, p2, c);
}

int H[4][4] = {
    {2, 1, -2, 1},
    {-3, -2, 3, -1},
    {0, 1, 0, 0},
    {1, 0, 0, 0}
};


void DrawHermite(HDC hdc, Point p[], int numOfSteps, COLORREF c) {
    Point m1 = (p[2] - p[0]) / 2, m2 = (p[3] - p[1]) / 2;

    Point input[] = {p[0], m1, p[3], m2};

    float step = 1.0 / numOfSteps;

    for (float i = 0; i <= 1; i += step) {
        Point *co = mul(H, input);
        Point res = co[0] * pow(i, 3) + co[1] * pow(i, 2) + co[2] * i + co[3];
        SetPixel(hdc, Round(res.x), Round(res.y), c);
    }
}

Point p1, p2, p3, p4;
int count = 0;

LRESULT WndProc(HWND hwnd, UINT m, WPARAM wp, LPARAM lp)
{
    HDC hdc;
    switch (m)
    {
        case WM_LBUTTONDOWN:
            count++;
            if (count == 1)
                p1 = Point(LOWORD(lp), HIWORD(lp));
            else if (count == 2)
                p2 = Point(LOWORD(lp), HIWORD(lp));
            else if (count == 3)
                p3 = Point(LOWORD(lp), HIWORD(lp));
            break;
        case WM_LBUTTONUP:
            if (count == 4) {
                p4 = Point(LOWORD(lp), HIWORD(lp));
                Point points[] = {p1, p2, p3, p4};
                hdc = GetDC(hwnd);
                DrawHermite(hdc, points, 3000, RGB(255, 0, 0));
                DrawBezier2(hdc, points, RGB(0, 255, 0));
                ReleaseDC(hwnd, hdc);
                count = 0;
            }
            break;
        case WM_CLOSE:
            DestroyWindow(hwnd); break;
        case WM_DESTROY:
            PostQuitMessage(0); break;
        default:return DefWindowProc(hwnd, m, wp, lp);
    }
    return 0;
}
int APIENTRY WinMain(HINSTANCE hi, HINSTANCE pi, LPSTR cmd, int nsh)
{
    WNDCLASS wc;
    wc.cbClsExtra = 0;
    wc.cbWndExtra = 0;
    wc.hbrBackground = (HBRUSH)GetStockObject(DKGRAY_BRUSH);
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hIcon = LoadIcon(NULL, IDI_WINLOGO);
    wc.lpszClassName = reinterpret_cast<LPCSTR>(L"MyClass");
    wc.lpszMenuName = NULL;
    wc.lpfnWndProc = WndProc;
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.hInstance = hi;
    RegisterClass(&wc);
    HWND hwnd = CreateWindow(reinterpret_cast<LPCSTR>(L"MyClass"), reinterpret_cast<LPCSTR>(L"Hello World!"), WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, NULL, NULL, hi, 0);
    ShowWindow(hwnd, nsh);
    UpdateWindow(hwnd);
    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0) > 0)
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return 0;
}
