#define UNICODE

#include <Windows.h>
#include <cmath>

using std::pair;

struct Point {
    int x, y;
};

union Outcode {
    struct {
        unsigned left: 1;
        unsigned right: 1;
        unsigned top: 1;
        unsigned bottom: 1;
    };
    unsigned all: 4;
};

void draw_line(const HDC hdc, const Point p1, const Point p2, const COLORREF c) {
    const int dx = p2.x - p1.x;
    const int dy = p2.y - p1.y;

    const int steps = abs(dx) > abs(dy) ? abs(dx) : abs(dy);

    const float x_inc = dx / static_cast<float>(steps);
    const float y_inc = dy / static_cast<float>(steps);

    float x = p1.x;
    float y = p1.y;

    for (int i = 0; i <= steps; ++i) {
        SetPixel(hdc, std::round(x), std::round(y), c);
        x += x_inc;
        y += y_inc;
    }
}

Outcode get_outcode(const Point p, const int xl, const int xr, const int yb, const int yt) {
    Outcode result {};

    result.all = 0;

    if (p.x < xl) {
        result.left = 1;
    }

    if (p.x > xr) {
        result.right = 1;
    }

    if (p.y < yb) {
        result.bottom = 1;
    }

    if (p.y > yt) {
        result.top = 1;
    }

    return result;
}

Point vIntersect(const Point p1, const Point p2, const int x_edge) {
    Point result {};

    result.x = x_edge;
    result.y = p1.y + (x_edge - p1.x) * (static_cast<double>(p2.y - p1.y)/(p2.x-p1.x));

    return result;
}

Point hIntersect(const Point p1, const Point p2, const int y_edge) {
    Point result {};

    result.y = y_edge;
    result.x = p1.x + (y_edge - p1.y) * (static_cast<double>(p2.x - p1.x)/(p2.y-p1.y));

    return result;
}

void clip_line(const HDC hdc, Point p1, Point p2, const int xl, const int xr, const int yb, const int yt) {
    Outcode out1 = get_outcode(p1, xl, xr, yb, yt);
    Outcode out2 = get_outcode(p2, xl, xr, yb, yt);

    while (true) {
        if (out1.all == 0 && out2.all == 0) {
            draw_line(hdc, p1, p2, RGB(0, 0, 0));
            return;
        }

        if ((out1.all & out2.all) != 0) {
            return;
        }

        Outcode out_code{};
        Point* point;

        if (out1.all != 0) {
            out_code = out1;
            point = &p1;
        } else {
            out_code = out2;
            point = &p2;
        }

        if (out_code.left) {
            *point = vIntersect(p1, p2, xl);
        }
        else if (out_code.right) {
            *point = vIntersect(p1, p2, xr);
        }
        else if (out_code.bottom) {
            *point = hIntersect(p1, p2, yb);
        }
        else if (out_code.top) {
            *point = hIntersect(p1, p2, yt);
        }

        if (point == &p1) {
            out1 = get_outcode(p1, xl, xr, yb, yt);
        } else {
            out2 = get_outcode(p2, xl, xr, yb, yt);
        }
    }
}

LRESULT WndProc(const HWND hwnd, const UINT m, const WPARAM wp, const LPARAM lp) {
    HDC hdc;
    static int left, top, right, bottom;
    static Point p1, p2;
    static bool draw = false;
    switch (m) {
    case WM_PAINT: {
        PAINTSTRUCT ps;
        hdc = BeginPaint(hwnd, &ps);

        RECT rect;
        GetClientRect(hwnd, &rect);
        const int width = rect.right - rect.left;
        const int height = rect.bottom - rect.top;

        constexpr int square_size {100};
        left = (width - square_size) / 2;
        top = (height - square_size) / 2;
        right = left + square_size;
        bottom = top + square_size;

        Rectangle(hdc, left, top, right + 1, bottom + 1);

        EndPaint(hwnd, &ps);
        break;
    }
    case WM_LBUTTONDOWN: {
        if (!draw) {
            p1 = {LOWORD(lp), HIWORD(lp)};
            draw = true;
        }
        else {
            hdc = GetDC(hwnd);
            p2 = {LOWORD(lp), HIWORD(lp)};
            clip_line(hdc, p1, p2, left, right, top, bottom);
            draw = false;
            ReleaseDC(hwnd, hdc);
        }
        break;
    }
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

int APIENTRY WinMain(const HINSTANCE hInstance, HINSTANCE, LPSTR, const int nShowCmd)
{
    WNDCLASS wc;
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
    const HWND hwnd = CreateWindow(L"MyClass", L"Hello World!", WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, nullptr, nullptr, hInstance, nullptr);
    ShowWindow(hwnd, nShowCmd);
    UpdateWindow(hwnd);
    MSG msg;
    while (GetMessage(&msg, nullptr, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return 0;
}
