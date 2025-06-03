#include <windows.h>
#include <commdlg.h>
#include <bits/stdc++.h>
#include <cmath>

using namespace std;

// Command IDs
enum MenuID {
    ID_COLOR_RED = 1,
    ID_COLOR_GREEN,
    ID_COLOR_BLUE,
    ID_CURSOR_ARROW,
    ID_CURSOR_CROSS,
    ID_CLEAR_SCREEN,
    ID_SAVE_IMAGE,
    ID_LOAD_IMAGE,
    ID_DDA_LINE,
    ID_MIDPOINT_LINE,
    ID_PARAMETRIC_LINE,
    ID_DIRECT_CIRCLE,
    ID_POLAR_CIRCLE,
    ID_ITER_POLAR_CIRCLE,
    ID_MIDPOINT_CIRCLE,
    ID_MOD_MIDPOINT_CIRCLE,
    ID_FILL_CIRCLE_LINES,
    ID_FILL_CIRCLE_CIRCLES,
    ID_FILL_SQUARE_HERMIT,
    ID_FILL_RECTANGLE_BEZIER,
    ID_FILL_CONVEX,
    ID_FILL_NONCONVEX,
    ID_FLOOD_RECURSIVE,
    ID_FLOOD_NONRECURSIVE,
    ID_CARDINAL_SPLINE,
    ID_ELLIPSE_DIRECT,
    ID_ELLIPSE_POLAR,
    ID_ELLIPSE_MIDPOINT,
    ID_CLIP_RECT_POINT,
    ID_CLIP_RECT_LINE,
    ID_CLIP_RECT_POLYGON,
    ID_CLIP_SQUARE_POINT,
    ID_CLIP_SQUARE_LINE
};

HBITMAP loadedBitmap = NULL;
COLORREF currentColor = RGB(0, 0, 0);
POINT points[4];
int pointCount = 0;
int currentAlgorithm = ID_DDA_LINE;
HCURSOR currentCursor = LoadCursor(NULL, IDC_ARROW);

// Utility to draw line using DDA
void DrawDDALine(HDC hdc, POINT p1, POINT p2, COLORREF color) {
    int dx = p2.x - p1.x;
    int dy = p2.y - p1.y;

    if (abs(dx) >= abs(dy)) {
        if (p1.x > p2.x) swap(p1, p2);
        float m = (float)dy / dx;
        float y = p1.y;
        for (int x = p1.x; x <= p2.x; ++x) {
            SetPixel(hdc, x, round(y), color);
            y += m;
        }
    } else {
        if (p1.y > p2.y) swap(p1, p2);
        float mInv = (float)dx / dy;
        float x = p1.x;
        for (int y = p1.y; y <= p2.y; ++y) {
            SetPixel(hdc, round(x), y, color);
            x += mInv;
        }
    }
}

// Midpoint Line Algorithm
void DrawMidpointLine(HDC hdc, POINT p1, POINT p2, COLORREF color) {
    int dx = p2.x - p1.x;
    int dy = p2.y - p1.y;
    int x = p1.x, y = p1.y;

    int xInc = (dx > 0) ? 1 : -1;
    int yInc = (dy > 0) ? 1 : -1;
    dx = abs(dx);
    dy = abs(dy);

    if (dx >= dy) { // Slope <= 1
        int d = 2 * dy - dx;
        int d1 = 2 * dy;
        int d2 = 2 * (dy - dx);

        SetPixel(hdc, x, y, color);
        while (x != p2.x) {
            if (d < 0) {
                d += d1;
            } else {
                d += d2;
                y += yInc;
            }
            x += xInc;
            SetPixel(hdc, x, y, color);
        }
    } else { // Slope > 1
        int d = 2 * dx - dy;
        int d1 = 2 * dx;
        int d2 = 2 * (dx - dy);

        SetPixel(hdc, x, y, color);
        while (y != p2.y) {
            if (d < 0) {
                d += d1;
            } else {
                d += d2;
                x += xInc;
            }
            y += yInc;
            SetPixel(hdc, x, y, color);
        }
    }
}

// Parametric Line Algorithm
void DrawParametricLine(HDC hdc, POINT p1, POINT p2, COLORREF color) {
    float t;
    for (t = 0; t <= 1; t += 0.001) {
        int x = p1.x + t * (p2.x - p1.x);
        int y = p1.y + t * (p2.y - p1.y);
        SetPixel(hdc, x, y, color);
    }
}

void Draw8Points(HDC hdc, POINT center, int x, int y, COLORREF color) {
    SetPixel(hdc, center.x + x, center.y + y, color);
    SetPixel(hdc, center.x - x, center.y + y, color);
    SetPixel(hdc, center.x + x, center.y - y, color);
    SetPixel(hdc, center.x - x, center.y - y, color);
    SetPixel(hdc, center.x + y, center.y + x, color);
    SetPixel(hdc, center.x - y, center.y + x, color);
    SetPixel(hdc, center.x + y, center.y - x, color);
    SetPixel(hdc, center.x - y, center.y - x, color);
}

// Direct Circle Algorithm (using Draw8Points)
void DrawDirectCircle(HDC hdc, POINT center, int radius, COLORREF color) {
    int x = 0, y = radius;
    int radiusSq = radius * radius;
    
    while (x <= y) {
        Draw8Points(hdc, center, x, y, color);
        x++;
        y = round(sqrt(radiusSq - x * x));
    }
}

// Polar Circle Algorithm (using Draw8Points)
void DrawPolarCircle(HDC hdc, POINT center, int radius, COLORREF color) {
    double dtheta = 1.0 / radius;
    for (double theta = 0; theta <= M_PI/4; theta += dtheta) {
        int x = round(radius * cos(theta));
        int y = round(radius * sin(theta));
        Draw8Points(hdc, center, x, y, color);
    }
}

// Iterative Polar Circle Algorithm (using Draw8Points)
void DrawIterativePolarCircle(HDC hdc, POINT center, int radius, COLORREF color) {
    double x = radius, y = 0;
    double dtheta = 1.0 / radius;
    double cdtheta = cos(dtheta), sdtheta = sin(dtheta);
    
    while (x >= y) {
        Draw8Points(hdc, center, round(x), round(y), color);
        
        // Rotate point
        double x_new = x * cdtheta - y * sdtheta;
        y = x * sdtheta + y * cdtheta;
        x = x_new;
    }
}

// Midpoint Circle Algorithm (using Draw8Points)
void DrawMidpointCircle(HDC hdc, POINT center, int radius, COLORREF color) {
    int x = 0, y = radius;
    int d = 1 - radius;
    
    while (x <= y) {
        Draw8Points(hdc, center, x, y, color);
        
        if (d < 0) {
            d += 2 * x + 3;
        } else {
            d += 2 * (x - y) + 5;
            y--;
        }
        x++;
    }
}

// Modified Midpoint Circle Algorithm (using Draw8Points)
void DrawModifiedMidpointCircle(HDC hdc, POINT center, int radius, COLORREF color) {
    int x = 0, y = radius;
    int d = 1 - radius;
    int d1 = 3;
    int d2 = -2 * radius + 5;
    
    while (x <= y) {
        Draw8Points(hdc, center, x, y, color);
        if (d < 0) {
            d += d1;
            d1 += 2;
            d2 += 2;
        } else {
            d += d2;
            d1 += 2;
            d2 += 4;
            y--;
        }
        x++;
    }
}

// Save screen content using BitBlt and file dialog
void SaveScreenToFile(HWND hwnd) {
    HDC hdcScreen = GetDC(hwnd);
    HDC hdcMem = CreateCompatibleDC(hdcScreen);
    RECT rc; 
    GetClientRect(hwnd, &rc);
    
    HBITMAP hBitmap = CreateCompatibleBitmap(hdcScreen, rc.right, rc.bottom);
    SelectObject(hdcMem, hBitmap);
    BitBlt(hdcMem, 0, 0, rc.right, rc.bottom, hdcScreen, 0, 0, SRCCOPY);

    OPENFILENAME ofn = { sizeof(ofn) };
    char fileName[MAX_PATH] = "untitled.bmp";  // Default name with extension
    
    ofn.lpstrFilter = "Bitmap Files (*.bmp)\0*.bmp\0All Files (*.*)\0*.*\0";
    ofn.lpstrFile = fileName;
    ofn.nMaxFile = MAX_PATH;
    ofn.lpstrDefExt = "bmp";  // Default extension if user doesn't specify
    ofn.lpstrTitle = "Save As Bitmap";
    ofn.Flags = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST;

    if (GetSaveFileName(&ofn)) {
        // Ensure the extension is .bmp
        char* ext = strrchr(fileName, '.');
        if (ext == NULL || _stricmp(ext, ".bmp") != 0) {
            if (strlen(fileName) + 4 < MAX_PATH) {  // Check buffer space
                strcat(fileName, ".bmp");
            } else {
                MessageBox(hwnd, "Filename too long", "Error", MB_OK | MB_ICONERROR);
                return;
            }
        }

        BITMAPFILEHEADER bmfHeader = {0};
        BITMAPINFOHEADER bi = {0};
        bi.biSize = sizeof(BITMAPINFOHEADER);
        bi.biWidth = rc.right;
        bi.biHeight = -rc.bottom;  // Negative for top-down DIB
        bi.biPlanes = 1;
        bi.biBitCount = 24;
        bi.biCompression = BI_RGB;

        DWORD dwBmpSize = ((rc.right * 3 + 3) & ~3) * abs(rc.bottom);  // 24-bit aligned row size
        BYTE* lpBits = new BYTE[dwBmpSize];

        GetDIBits(hdcMem, hBitmap, 0, rc.bottom, lpBits, (BITMAPINFO*)&bi, DIB_RGB_COLORS);

        ofstream file(fileName, ios::binary);
        if (file.is_open()) {
            bmfHeader.bfType = 0x4D42;  // 'BM'
            bmfHeader.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
            bmfHeader.bfSize = bmfHeader.bfOffBits + dwBmpSize;
            
            file.write((char*)&bmfHeader, sizeof(bmfHeader));
            file.write((char*)&bi, sizeof(bi));
            file.write((char*)lpBits, dwBmpSize);
            file.close();
        } else {
            MessageBox(hwnd, "Failed to create file", "Error", MB_OK | MB_ICONERROR);
        }

        delete[] lpBits;
    }

    DeleteObject(hBitmap);
    DeleteDC(hdcMem);
    ReleaseDC(hwnd, hdcScreen);
}

void LoadImageFromFile(HWND hwnd) {
    OPENFILENAME ofn = { sizeof(ofn) };
    char fileName[MAX_PATH] = "";
    
    ofn.hwndOwner = hwnd;
    ofn.lpstrFile = fileName;
    ofn.nMaxFile = MAX_PATH;
    ofn.lpstrFilter = "Bitmap Files (*.bmp)\0*.bmp\0All Files (*.*)\0*.*\0";
    ofn.lpstrTitle = "Select Bitmap File to Open";
    ofn.Flags = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST | OFN_HIDEREADONLY;
    
    if (GetOpenFileName(&ofn)) {
        // Try to open the file to verify it exists
        HANDLE hFile = CreateFile(
            fileName,
            GENERIC_READ,
            FILE_SHARE_READ,
            NULL,
            OPEN_EXISTING,
            FILE_ATTRIBUTE_NORMAL,
            NULL
        );
        
        if (hFile == INVALID_HANDLE_VALUE) {
            DWORD err = GetLastError();
            char msg[256];
            sprintf(msg, "Cannot open file (Error %d)", err);
            MessageBox(hwnd, msg, "Error", MB_OK | MB_ICONERROR);
            return;
        }
        CloseHandle(hFile);

        // Load the bitmap
        HBITMAP hBitmap = (HBITMAP)LoadImage(
            NULL,
            fileName,
            IMAGE_BITMAP,
            0,
            0,
            LR_LOADFROMFILE | LR_CREATEDIBSECTION | LR_DEFAULTSIZE
        );

        if (!hBitmap) {
            DWORD err = GetLastError();
            char msg[256];
            sprintf(msg, "Failed to load image (Error %d)", err);
            MessageBox(hwnd, msg, "Error", MB_OK | MB_ICONERROR);
            return;
        }

        // Clean up previous bitmap if exists
        if (loadedBitmap) {
            DeleteObject(loadedBitmap);
        }

        loadedBitmap = hBitmap;
        InvalidateRect(hwnd, NULL, TRUE);
    }
}

void AddMenus(HWND hwnd) {
    HMENU hMain = CreateMenu();
    HMENU hColor = CreateMenu(), hCursor = CreateMenu(), hLines = CreateMenu(), hCircles = CreateMenu();
    HMENU hFills = CreateMenu(), hClipping = CreateMenu(), hEllipse = CreateMenu();

    // Color menu
    AppendMenu(hColor, MF_STRING, ID_COLOR_RED, "Red");
    AppendMenu(hColor, MF_STRING, ID_COLOR_GREEN, "Green");
    AppendMenu(hColor, MF_STRING, ID_COLOR_BLUE, "Blue");

    // Cursor menu
    AppendMenu(hCursor, MF_STRING, ID_CURSOR_ARROW, "Arrow");
    AppendMenu(hCursor, MF_STRING, ID_CURSOR_CROSS, "Cross");

    // Lines menu
    AppendMenu(hLines, MF_STRING, ID_DDA_LINE, "DDA Line");
    AppendMenu(hLines, MF_STRING, ID_MIDPOINT_LINE, "Midpoint Line");
    AppendMenu(hLines, MF_STRING, ID_PARAMETRIC_LINE, "Parametric Line");

    // Other menus remain the same...
    AppendMenu(hCircles, MF_STRING, ID_DIRECT_CIRCLE, "Direct Circle");
    AppendMenu(hCircles, MF_STRING, ID_POLAR_CIRCLE, "Polar Circle");
    AppendMenu(hCircles, MF_STRING, ID_ITER_POLAR_CIRCLE, "Iterative Polar");
    AppendMenu(hCircles, MF_STRING, ID_MIDPOINT_CIRCLE, "Midpoint Circle");
    AppendMenu(hCircles, MF_STRING, ID_MOD_MIDPOINT_CIRCLE, "Modified Midpoint");

    AppendMenu(hFills, MF_STRING, ID_FILL_CIRCLE_LINES, "Fill Circle w/ Lines");
    AppendMenu(hFills, MF_STRING, ID_FILL_CIRCLE_CIRCLES, "Fill Circle w/ Circles");
    AppendMenu(hFills, MF_STRING, ID_FILL_SQUARE_HERMIT, "Fill Square Hermite");
    AppendMenu(hFills, MF_STRING, ID_FILL_RECTANGLE_BEZIER, "Fill Rect Bezier");
    AppendMenu(hFills, MF_STRING, ID_FILL_CONVEX, "Fill Convex");
    AppendMenu(hFills, MF_STRING, ID_FILL_NONCONVEX, "Fill Non-Convex");
    AppendMenu(hFills, MF_STRING, ID_FLOOD_RECURSIVE, "Flood Fill Recursive");
    AppendMenu(hFills, MF_STRING, ID_FLOOD_NONRECURSIVE, "Flood Fill Non-Recursive");
    AppendMenu(hFills, MF_STRING, ID_CARDINAL_SPLINE, "Cardinal Spline");

    AppendMenu(hEllipse, MF_STRING, ID_ELLIPSE_DIRECT, "Direct Ellipse");
    AppendMenu(hEllipse, MF_STRING, ID_ELLIPSE_POLAR, "Polar Ellipse");
    AppendMenu(hEllipse, MF_STRING, ID_ELLIPSE_MIDPOINT, "Midpoint Ellipse");

    AppendMenu(hClipping, MF_STRING, ID_CLIP_RECT_POINT, "Clip Point (Rect)");
    AppendMenu(hClipping, MF_STRING, ID_CLIP_RECT_LINE, "Clip Line (Rect)");
    AppendMenu(hClipping, MF_STRING, ID_CLIP_RECT_POLYGON, "Clip Polygon (Rect)");
    AppendMenu(hClipping, MF_STRING, ID_CLIP_SQUARE_POINT, "Clip Point (Square)");
    AppendMenu(hClipping, MF_STRING, ID_CLIP_SQUARE_LINE, "Clip Line (Square)");

    // Main menu organization
    AppendMenu(hMain, MF_POPUP, (UINT_PTR)hColor, "Color");
    AppendMenu(hMain, MF_POPUP, (UINT_PTR)hCursor, "Cursor");
    AppendMenu(hMain, MF_POPUP, (UINT_PTR)hLines, "Line");
    AppendMenu(hMain, MF_POPUP, (UINT_PTR)hCircles, "Circle");
    AppendMenu(hMain, MF_POPUP, (UINT_PTR)hFills, "Fills");
    AppendMenu(hMain, MF_POPUP, (UINT_PTR)hEllipse, "Ellipse");
    AppendMenu(hMain, MF_POPUP, (UINT_PTR)hClipping, "Clipping");

    AppendMenu(hMain, MF_STRING, ID_CLEAR_SCREEN, "Clear");
    AppendMenu(hMain, MF_STRING, ID_SAVE_IMAGE, "Save Screen");
    AppendMenu(hMain, MF_STRING, ID_LOAD_IMAGE, "Load Image");

    SetMenu(hwnd, hMain);
}

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    switch (msg) {
    case WM_COMMAND:
        switch (LOWORD(wParam)) {
        case ID_COLOR_RED:   currentColor = RGB(255, 0, 0); break;
        case ID_COLOR_GREEN: currentColor = RGB(0, 255, 0); break;
        case ID_COLOR_BLUE:  currentColor = RGB(0, 0, 255); break;

        case ID_CURSOR_ARROW:
            currentCursor = LoadCursor(NULL, IDC_ARROW);
            SetCursor(currentCursor);
            break;

        case ID_CURSOR_CROSS:
            currentCursor = LoadCursor(NULL, IDC_CROSS);
            SetCursor(currentCursor);
            break;

        case ID_DDA_LINE:
        case ID_MIDPOINT_LINE:
        case ID_PARAMETRIC_LINE:
        case ID_DIRECT_CIRCLE:
        case ID_POLAR_CIRCLE:
        case ID_ITER_POLAR_CIRCLE:
        case ID_MIDPOINT_CIRCLE:
        case ID_MOD_MIDPOINT_CIRCLE:
        case ID_FILL_CIRCLE_LINES:
        case ID_FILL_CIRCLE_CIRCLES:
        case ID_FILL_SQUARE_HERMIT:
        case ID_FILL_RECTANGLE_BEZIER:
        case ID_FILL_CONVEX:
        case ID_FILL_NONCONVEX:
        case ID_FLOOD_RECURSIVE:
        case ID_FLOOD_NONRECURSIVE:
        case ID_CARDINAL_SPLINE:
        case ID_ELLIPSE_DIRECT:
        case ID_ELLIPSE_POLAR:
        case ID_ELLIPSE_MIDPOINT:
        case ID_CLIP_RECT_POINT:
        case ID_CLIP_RECT_LINE:
        case ID_CLIP_RECT_POLYGON:
        case ID_CLIP_SQUARE_POINT:
        case ID_CLIP_SQUARE_LINE:
            // Clear points when switching algorithms
            pointCount = 0;
            currentAlgorithm = LOWORD(wParam);
            cout << "Algorithm changed to: " << currentAlgorithm << endl;
            break;

        case ID_CLEAR_SCREEN:
            pointCount = 0;
            if (loadedBitmap) {
                DeleteObject(loadedBitmap);
                loadedBitmap = NULL;
            }
            InvalidateRect(hwnd, NULL, TRUE);
            break;

        case ID_SAVE_IMAGE:
            pointCount = 0;
            SaveScreenToFile(hwnd);
            break;

        case ID_LOAD_IMAGE:
            pointCount = 0;
            LoadImageFromFile(hwnd);
            break;
        }
        break;

    case WM_SETCURSOR:
        if (LOWORD(lParam) == HTCLIENT) {
            SetCursor(currentCursor);
            return TRUE;
        }
        break;

    case WM_LBUTTONDOWN: {
        POINT pt = { LOWORD(lParam), HIWORD(lParam) };
        
        if (pointCount < 4) {
            points[pointCount++] = pt;
        }

        // For line algorithms, draw when we have 2 points
        if (pointCount == 2 && 
            (currentAlgorithm == ID_DDA_LINE || 
             currentAlgorithm == ID_MIDPOINT_LINE || 
             currentAlgorithm == ID_PARAMETRIC_LINE)) {
            
            HDC hdc = GetDC(hwnd);
            switch (currentAlgorithm) {
                case ID_DDA_LINE:
                    DrawDDALine(hdc, points[0], points[1], currentColor);
                    break;
                case ID_MIDPOINT_LINE:
                    DrawMidpointLine(hdc, points[0], points[1], currentColor);
                    break;
                case ID_PARAMETRIC_LINE:
                    DrawParametricLine(hdc, points[0], points[1], currentColor);
                    break;
            }
            ReleaseDC(hwnd, hdc);
            pointCount = 0; // Reset points for next shape
        } else if (pointCount == 2 && 
            (currentAlgorithm == ID_DIRECT_CIRCLE ||
            currentAlgorithm == ID_POLAR_CIRCLE ||
            currentAlgorithm == ID_ITER_POLAR_CIRCLE ||
            currentAlgorithm == ID_MIDPOINT_CIRCLE ||
            currentAlgorithm == ID_MOD_MIDPOINT_CIRCLE)) {
                HDC hdc = GetDC(hwnd);
                int radius = static_cast<int>(sqrt(
                    pow(points[1].x - points[0].x, 2) + 
                    pow(points[1].y - points[0].y, 2)
                ));
                
                switch (currentAlgorithm) {
                    case ID_DIRECT_CIRCLE:
                        DrawDirectCircle(hdc, points[0], radius, currentColor);
                        break;
                    case ID_POLAR_CIRCLE:
                        DrawPolarCircle(hdc, points[0], radius, currentColor);
                        break;
                    case ID_ITER_POLAR_CIRCLE:
                        DrawIterativePolarCircle(hdc, points[0], radius, currentColor);
                        break;
                    case ID_MIDPOINT_CIRCLE:
                        DrawMidpointCircle(hdc, points[0], radius, currentColor);
                        break;
                    case ID_MOD_MIDPOINT_CIRCLE:
                        DrawModifiedMidpointCircle(hdc, points[0], radius, currentColor);
                        break;
                }
                ReleaseDC(hwnd, hdc);
                pointCount = 0;
    }

        //TODO: Add the other algorithms.
        
        break;
    }

    case WM_PAINT: {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hwnd, &ps);

        if (loadedBitmap) {
            HDC hdcMem = CreateCompatibleDC(hdc);
            HBITMAP oldBitmap = (HBITMAP)SelectObject(hdcMem, loadedBitmap);
            BITMAP bmp;
            GetObject(loadedBitmap, sizeof(BITMAP), &bmp);
            BitBlt(hdc, 0, 0, bmp.bmWidth, bmp.bmHeight, hdcMem, 0, 0, SRCCOPY);
            SelectObject(hdcMem, oldBitmap);
            DeleteDC(hdcMem);
        } 
        else {
            RECT rect;
            GetClientRect(hwnd, &rect);
            HBRUSH whiteBrush = (HBRUSH)GetStockObject(WHITE_BRUSH);
            FillRect(hdc, &rect, whiteBrush);
        }

        EndPaint(hwnd, &ps);
        break;
    }

    case WM_DESTROY:
        if (loadedBitmap) DeleteObject(loadedBitmap);
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hwnd, msg, wParam, lParam);
    }
    return 0;
}

int WINAPI WinMain(HINSTANCE hInst, HINSTANCE, LPSTR, int nCmdShow) {
    AllocConsole(); freopen("CONOUT$", "w", stdout);
    cout << "Program started.\n";

    WNDCLASS wc = { 0 };
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInst;
    wc.lpszClassName = "GraphicsApp";
    wc.hCursor = currentCursor;
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);

    RegisterClass(&wc);

    HWND hwnd = CreateWindow("GraphicsApp", "Graphics Application", WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, 1000, 700, NULL, NULL, hInst, NULL);
    AddMenus(hwnd);
    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);

    MSG msg;
    while (GetMessage(&msg, 0, 0, 0)) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return 0;
}
