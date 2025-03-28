#include <windows.h>
#include <wingdi.h>
#include <math.h>
#ifdef UNICODE
#define WinMain wWinMain
#endif

int xc, yc;

void DrawPoints(HDC hdc, int xc, int yc, int x, int y, COLORREF c)
{
    SetPixel(hdc, xc + x, yc + y, c);
    SetPixel(hdc, xc - x, yc + y, c);
    SetPixel(hdc, xc + x, yc - y, c);
    SetPixel(hdc, xc - x, yc - y, c);
    SetPixel(hdc, xc + y, yc + x, c);
    SetPixel(hdc, xc - y, yc + x, c);
    SetPixel(hdc, xc - y, yc - x, c);
    SetPixel(hdc, xc + y, yc - x, c);
}


void DrawCircleBres(HDC hdc, int xc, int yc, int R, COLORREF c)
{
    int x = 0, y = R;
    int d = 1 - R;
    int d1 = 3;
    int d2 = -2 * R + 5;


    DrawPoints(hdc, xc, yc, x, y, c);
    while (x <= y)
    {
        x++;
        if (d <= 0)
        {
            d += d1;
            d2 += 2;
        }
        else
        {
            d += d2;
            d2 += 4;
            y--;
        }
        d1 += 2;
        DrawPoints(hdc, xc, yc, x, y, c);
    }
}
LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp)
{
    HDC hdc;
    int x, y;
    switch (msg)
    {
        case WM_LBUTTONDOWN:
            xc = LOWORD(lp);
            yc = HIWORD(lp);
            break;
        case WM_LBUTTONUP:
            hdc = GetDC(hwnd);
            x = LOWORD(lp);
            y = HIWORD(lp);
            int R = (int)sqrt(pow((xc-x),2) + pow((yc-y),2));
            DrawCircleBres(hdc, xc, yc, R, RGB(255, 0, 0));
            ReleaseDC(hwnd, hdc);
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

int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
    WNDCLASSW wc = {0};
    wc.hbrBackground = CreateSolidBrush(RGB(50, 50, 50));
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
    wc.lpszClassName = L"MyClass";
    wc.lpfnWndProc = WndProc;
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.hInstance = hInstance;
    if (!RegisterClassW(&wc))
    {
        MessageBoxW(NULL, L"Window Registration Failed!", L"Error", MB_ICONERROR);
        return 0;
    }
    HWND hwnd = CreateWindowW(L"MyClass", L"Circle Drawing - Polar Algorithm", WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, 500, 500, NULL, NULL, hInstance, NULL);
    if (!hwnd)
    {
        MessageBoxW(NULL, L"Window Creation Failed!", L"Error", MB_ICONERROR);
        return 0;
    }
    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);
    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0))
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return (int)msg.wParam;
}
