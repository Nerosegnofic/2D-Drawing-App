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
COLORREF currentColor = RGB(255, 0, 0);
COLORREF bgColor = RGB(255, 255, 255);
POINT points[200];
int pointCount = 0;
int currentAlgorithm = ID_DDA_LINE;
HCURSOR currentCursor = LoadCursor(NULL, IDC_ARROW);

typedef struct {
    int xleft, xright;
} EdgeTable[80000];

struct Node {
    double x, minv;
    int ymax;
    Node(double _x = 0, int _ymax = 0, double _minv = 0) : x(_x), minv(_minv), ymax(_ymax) {}

    bool operator<(const Node& other) const {
        return x < other.x;
    }
};

struct square_clip_window {
    int left, right, top, bottom;
    int size;
};

struct rect_clip_window {
    int left, right, top, bottom;
    int width, height;
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

square_clip_window squareClipWindow = {0, 0, 0, 0, 100}; // Default 100x100 square
rect_clip_window rectClipWindow = {0, 0, 0, 0, 150, 100}; // Default 150x100 rectangle

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

// Initialize clipping window in the center of the screen
void InitializeSquareClipWindow(HWND hwnd) {
    RECT rect;
    GetClientRect(hwnd, &rect);
    int width = rect.right - rect.left;
    int height = rect.bottom - rect.top;

    squareClipWindow.left = (width - squareClipWindow.size) / 2;
    squareClipWindow.top = (height - squareClipWindow.size) / 2;
    squareClipWindow.right = squareClipWindow.left + squareClipWindow.size;
    squareClipWindow.bottom = squareClipWindow.top + squareClipWindow.size;
}

void InitializeRectClipWindow(HWND hwnd) {
    RECT rect;
    GetClientRect(hwnd, &rect);
    int screenWidth = rect.right - rect.left;
    int screenHeight = rect.bottom - rect.top;

    rectClipWindow.left = (screenWidth - rectClipWindow.width) / 2;
    rectClipWindow.top = (screenHeight - rectClipWindow.height) / 2;
    rectClipWindow.right = rectClipWindow.left + rectClipWindow.width;
    rectClipWindow.bottom = rectClipWindow.top + rectClipWindow.height;
}

// Draw the clipping window
void DrawSquareClipWindow(HDC hdc) {
    HPEN pen = CreatePen(PS_SOLID, 2, RGB(255, 0, 0)); // Red border
    HPEN oldPen = (HPEN)SelectObject(hdc, pen);

    Rectangle(hdc, squareClipWindow.left, squareClipWindow.top, squareClipWindow.right + 1, squareClipWindow.bottom + 1);

    SelectObject(hdc, oldPen);
    DeleteObject(pen);
}

void DrawRectClipWindow(HDC hdc) {
    HPEN pen = CreatePen(PS_SOLID, 2, RGB(0, 0, 255)); // Blue border to distinguish from square
    HPEN oldPen = (HPEN)SelectObject(hdc, pen);

    Rectangle(hdc, rectClipWindow.left, rectClipWindow.top,
              rectClipWindow.right + 1, rectClipWindow.bottom + 1);

    SelectObject(hdc, oldPen);
    DeleteObject(pen);
}

// Generate outcode for a point (from your algorithm)
Outcode GetOutcode(POINT p, int xl, int xr, int yb, int yt) {
    Outcode result = {0};

    if (p.x < xl) result.left = 1;
    if (p.x > xr) result.right = 1;
    if (p.y < yb) result.bottom = 1;
    if (p.y > yt) result.top = 1;

    return result;
}

// Vertical intersection (from your algorithm)
POINT VIntersect(POINT p1, POINT p2, int x_edge) {
    POINT result;
    result.x = x_edge;
    result.y = p1.y + (x_edge - p1.x) * (double)(p2.y - p1.y) / (p2.x - p1.x);
    return result;
}

// Horizontal intersection (from your algorithm)
POINT HIntersect(POINT p1, POINT p2, int y_edge) {
    POINT result;
    result.y = y_edge;
    result.x = p1.x + (y_edge - p1.y) * (double)(p2.x - p1.x) / (p2.y - p1.y);
    return result;
}

// Point clipping using square window
bool ClipPointSquare(HDC hdc, POINT p, COLORREF color) {
    // Check if point is inside the clipping window
    if (p.x >= squareClipWindow.left && p.x <= squareClipWindow.right &&
        p.y >= squareClipWindow.top && p.y <= squareClipWindow.bottom) {
        // Point is inside, draw it
        SetPixel(hdc, p.x, p.y, color);
        return true;
    }
    return false; // Point is clipped (outside window)
}

bool ClipPointRect(HDC hdc, POINT p, COLORREF color) {
    // Check if point is inside the rectangle clipping window
    if (p.x >= rectClipWindow.left && p.x <= rectClipWindow.right &&
        p.y >= rectClipWindow.top && p.y <= rectClipWindow.bottom) {
        // Point is inside, draw it
        SetPixel(hdc, p.x, p.y, color);
        return true;
        }
    return false; // Point is clipped (outside window)
}

// Line clipping using Cohen-Sutherland algorithm with square window (adapted from your code)
void ClipLineSquare(HDC hdc, POINT p1, POINT p2, COLORREF color) {
    Outcode out1 = GetOutcode(p1, squareClipWindow.left, squareClipWindow.right, squareClipWindow.top, squareClipWindow.bottom);
    Outcode out2 = GetOutcode(p2, squareClipWindow.left, squareClipWindow.right, squareClipWindow.top, squareClipWindow.bottom);

    while (true) {
        // Both points inside window
        if (out1.all == 0 && out2.all == 0) {
            DrawDDALine(hdc, p1, p2, color);
            return;
        }

        // Both points outside window on same side
        if ((out1.all & out2.all) != 0) {
            return; // Line is completely outside
        }

        // Line crosses window boundary
        Outcode outcode;
        POINT* point;

        if (out1.all != 0) {
            outcode = out1;
            point = &p1;
        } else {
            outcode = out2;
            point = &p2;
        }

        // Find intersection point
        if (outcode.left) {
            *point = VIntersect(p1, p2, squareClipWindow.left);
        } else if (outcode.right) {
            *point = VIntersect(p1, p2, squareClipWindow.right);
        } else if (outcode.bottom) {
            *point = HIntersect(p1, p2, squareClipWindow.top);
        } else if (outcode.top) {
            *point = HIntersect(p1, p2, squareClipWindow.bottom);
        }

        // Update outcode for the modified point
        if (point == &p1) {
            out1 = GetOutcode(p1, squareClipWindow.left, squareClipWindow.right, squareClipWindow.top, squareClipWindow.bottom);
        } else {
            out2 = GetOutcode(p2, squareClipWindow.left, squareClipWindow.right, squareClipWindow.top, squareClipWindow.bottom);
        }
    }
}

void ClipLineRect(HDC hdc, POINT p1, POINT p2, COLORREF color) {
    Outcode out1 = GetOutcode(p1, rectClipWindow.left, rectClipWindow.right,
                              rectClipWindow.top, rectClipWindow.bottom);
    Outcode out2 = GetOutcode(p2, rectClipWindow.left, rectClipWindow.right,
                              rectClipWindow.top, rectClipWindow.bottom);

    while (true) {
        // Both points inside window
        if (out1.all == 0 && out2.all == 0) {
            DrawDDALine(hdc, p1, p2, color);
            return;
        }

        // Both points outside window on same side
        if ((out1.all & out2.all) != 0) {
            return; // Line is completely outside
        }

        // Line crosses window boundary
        Outcode outcode;
        POINT* point;

        if (out1.all != 0) {
            outcode = out1;
            point = &p1;
        } else {
            outcode = out2;
            point = &p2;
        }

        // Find intersection point
        if (outcode.left) {
            *point = VIntersect(p1, p2, rectClipWindow.left);
        } else if (outcode.right) {
            *point = VIntersect(p1, p2, rectClipWindow.right);
        } else if (outcode.bottom) {
            *point = HIntersect(p1, p2, rectClipWindow.top);
        } else if (outcode.top) {
            *point = HIntersect(p1, p2, rectClipWindow.bottom);
        }

        // Update outcode for the modified point
        if (point == &p1) {
            out1 = GetOutcode(p1, rectClipWindow.left, rectClipWindow.right,
                              rectClipWindow.top, rectClipWindow.bottom);
        } else {
            out2 = GetOutcode(p2, rectClipWindow.left, rectClipWindow.right,
                              rectClipWindow.top, rectClipWindow.bottom);
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

void DrawDirectEllipse(HDC hdc, POINT center, int rx, int ry, COLORREF color) {
    if (rx <= 0 || ry <= 0) return;

    double rx2 = double(rx) * rx;
    double ry2 = double(ry) * ry;

    for (int xi = 0; xi <= rx; xi++) {
        double xd = double(xi);
        double yd = ry * sqrt(1.0 - (xd*xd) / rx2);
        int yi = int(round(yd));

        SetPixel(hdc, center.x + xi, center.y + yi, color);
        SetPixel(hdc, center.x - xi, center.y + yi, color);
        SetPixel(hdc, center.x + xi, center.y - yi, color);
        SetPixel(hdc, center.x - xi, center.y - yi, color);
    }

    for (int yi = 0; yi <= ry; yi++) {
        double yd = double(yi);
        double xd = rx * sqrt(1.0 - (yd*yd) / ry2);
        int xi = int(round(xd));

        SetPixel(hdc, center.x + xi, center.y + yi, color);
        SetPixel(hdc, center.x - xi, center.y + yi, color);
        SetPixel(hdc, center.x + xi, center.y - yi, color);
        SetPixel(hdc, center.x - xi, center.y - yi, color);
    }
}

// Polar Ellipse Algorithm
void DrawPolarEllipse(HDC hdc, POINT center, int rx, int ry, COLORREF color) {
    if (rx <= 0 || ry <= 0) return;  // Safety check

    double dtheta = 1.0 / max(rx, ry);
    for (double theta = 0; theta <= M_PI/2; theta += dtheta) {
        int x = round(rx * cos(theta));
        int y = round(ry * sin(theta));
        SetPixel(hdc, center.x + x, center.y + y, color);
        SetPixel(hdc, center.x - x, center.y + y, color);
        SetPixel(hdc, center.x + x, center.y - y, color);
        SetPixel(hdc, center.x - x, center.y - y, color);
    }
}

// Midpoint Ellipse Algorithm
void DrawMidpointEllipse(HDC hdc, POINT center, int rx, int ry, COLORREF color) {
    if (rx <= 0 || ry <= 0) return;  // Safety check

    double rx2 = rx * rx;
    double ry2 = ry * ry;
    double x = 0, y = ry;
    double p = ry2 - rx2 * ry + 0.25 * rx2;

    // Region 1
    while (ry2 * x < rx2 * y) {
        SetPixel(hdc, center.x + x, center.y + y, color);
        SetPixel(hdc, center.x - x, center.y + y, color);
        SetPixel(hdc, center.x + x, center.y - y, color);
        SetPixel(hdc, center.x - x, center.y - y, color);

        if (p < 0) {
            p += 2 * ry2 * x + ry2;
        } else {
            y--;
            p += 2 * ry2 * x - 2 * rx2 * y + ry2;
        }
        x++;
    }

    // Region 2
    p = ry2 * (x + 0.5) * (x + 0.5) + rx2 * (y - 1) * (y - 1) - rx2 * ry2;
    while (y >= 0) {
        SetPixel(hdc, center.x + x, center.y + y, color);
        SetPixel(hdc, center.x - x, center.y + y, color);
        SetPixel(hdc, center.x + x, center.y - y, color);
        SetPixel(hdc, center.x - x, center.y - y, color);

        if (p > 0) {
            p += -2 * rx2 * y + rx2;
        } else {
            x++;
            p += 2 * ry2 * x - 2 * rx2 * y + rx2;
        }
        y--;
    }
}

void FloodFillRecursive(HDC hdc, int x, int y, COLORREF targetColor, COLORREF fillColor) {
    // Get current pixel color
    COLORREF currentColor = GetPixel(hdc, x, y);

    // Base case: if current color is not target color or already filled
    if (currentColor != targetColor || currentColor == fillColor) {
        return;
    }

    // Fill current pixel
    SetPixel(hdc, x, y, fillColor);

    // Recursively fill in 4 directions (4-connected)
    FloodFillRecursive(hdc, x + 1, y, targetColor, fillColor);
    FloodFillRecursive(hdc, x - 1, y, targetColor, fillColor);
    FloodFillRecursive(hdc, x, y + 1, targetColor, fillColor);
    FloodFillRecursive(hdc, x, y - 1, targetColor, fillColor);
}

void FloodFillNonRecursive(HDC hdc, int startX, int startY, COLORREF targetColor, COLORREF fillColor) {
    // Use stack to simulate recursion
    stack<POINT> pixelStack;

    // Push starting point
    POINT startPoint = {startX, startY};
    pixelStack.push(startPoint);

    while (!pixelStack.empty()) {
        POINT currentPoint = pixelStack.top();
        pixelStack.pop();

        int x = currentPoint.x;
        int y = currentPoint.y;


        // Get current pixel color
        COLORREF currentColor = GetPixel(hdc, x, y);

        // Skip if not target color or already filled
        if (currentColor != targetColor || currentColor == fillColor) {
            continue;
        }

        // Fill current pixel
        SetPixel(hdc, x, y, fillColor);

        // Add neighboring pixels to stack (4-connected)
        POINT neighbors[4] = {
            {x + 1, y},     // Right
            {x - 1, y},     // Left
            {x, y + 1},     // Down
            {x, y - 1}      // Up
        };

        for (int i = 0; i < 4; i++) {
            pixelStack.push(neighbors[i]);
        }
    }
}

void DrawHermiteVerticalCurvePoints(vector<int>& pixelsY, int y0, int y1, int t0y, int t1y) {
    double a0 = y0;
    double a1 = t0y;
    double a2 = -3 * y0 + 3 * y1 - 2 * t0y - t1y;
    double a3 = 2 * y0 - 2 * y1 + t0y + t1y;
    pixelsY.clear();
    for (double t = 0; t <= 1.0; t += 0.001) {
        double y = a0 + a1 * t + a2 * t * t + a3 * t * t * t;
        pixelsY.push_back((int)y);
    }
}

void FillSquareWithVerticalHermite(HDC hdc, POINT p1, POINT p2, COLORREF color) {
    int width = abs(p2.x - p1.x);
    int height = abs(p2.y - p1.y);
    int side = min(width, height);
    int dirX = (p2.x > p1.x) ? 1 : -1;
    int dirY = (p2.y > p1.y) ? 1 : -1;
    int x0 = p1.x;
    int y0 = p1.y;
    int x1 = x0 + dirX * side;
    int y1 = y0 + dirY * side;
    int t0y = -side / 4;
    int t1y = -side / 4;
    int startX = min(x0, x1);
    int endX = max(x0, x1);
    int startY = min(y0, y1);
    int endY = max(y0, y1);
    int totalSteps = 1000;
    vector<vector<int>> allY(endX - startX + 1);
    for (int xi = 0; xi <= endX - startX; xi++) {
        int x = startX + xi;
        vector<int>& pixelsY = allY[xi];
        DrawHermiteVerticalCurvePoints(pixelsY, endY, startY, t0y, t1y);
    }
    for (int step = 0; step < totalSteps; step++) {
        for (int xi = 0; xi <= endX - startX; xi++) {
            int x = startX + xi;
            const vector<int>& pixelsY = allY[xi];
            if (step < pixelsY.size()) {
                SetPixel(hdc, x, pixelsY[step], color);
            }
        }
    }
}

void DrawBezierHorizontalCurve(HDC hdc, int x, POINT p0, POINT p1, POINT p2, POINT p3, COLORREF color) {
    double t, y;

    MoveToEx(hdc, x, p0.y, NULL);
    for (t = 0.0; t <= 1.0; t += 0.01) {
        y = pow(1 - t, 3) * p0.y +
            3 * pow(1 - t, 2) * t * p1.y +
            3 * (1 - t) * t * t * p2.y +
            t * t * t * p3.y;

        LineTo(hdc, x, (int)y);
    }
}

void FillRectangleWithHorizontalBezier(HDC hdc, POINT p1, POINT p2, COLORREF color) {
    int width = abs(p2.x - p1.x);
    int height = abs(p2.y - p1.y);

    int dirX = (p2.x > p1.x) ? 1 : -1;
    int dirY = (p2.y > p1.y) ? 1 : -1;

    int xStart = p1.x;
    int xEnd = p1.x + dirX * width;

    POINT top    = { 0, p1.y };
    POINT ctrl1  = { 0, p1.y + dirY * height / 3 };
    POINT ctrl2  = { 0, p1.y + dirY * 2 * height / 3 };
    POINT bottom = { 0, p1.y + dirY * height };

    HPEN pen = CreatePen(PS_SOLID, 1, color);
    HPEN oldPen = (HPEN)SelectObject(hdc, pen);

    for (int x = min(xStart, xEnd); x <= max(xStart, xEnd); x++) {
        top.x = ctrl1.x = ctrl2.x = bottom.x = x;
        DrawBezierHorizontalCurve(hdc, x, top, ctrl1, ctrl2, bottom, color);
    }

    SelectObject(hdc, oldPen);
    DeleteObject(pen);
}

void init(EdgeTable tbl) {
    for (int i = 0; i < 80000; i++) {
        tbl[i].xleft = 10000;
        tbl[i].xright = -10000;
    }
}

void edge2table(POINT &p1, POINT &p2, EdgeTable tbl){
    if (p1.y == p2.y) return;
    if (p1.y > p2.y) swap(p1, p2);
    int y = p1.y;
    double x = p1.x;
    double d = (double) (p2.x - p1.x) / (p2.y - p1.y);
    while (y < p2.y) {
        if (x < tbl[y].xleft) tbl[y].xleft = (int) ceil(x);
        if (x > tbl[y].xright) tbl[y].xright = (int) floor(x);
        y++;
        x += d;
    }
}

void polygon2table(POINT p[], int n, EdgeTable tbl) {
    POINT point = p[n - 1];
    for (int i = 0; i < n; i++) {
        POINT nextPoint = p[i];
        edge2table(point, nextPoint, tbl);

        point = p[i];
    }
}

void table2screen(HDC hdc, EdgeTable tbl, COLORREF color) {
    for (int y = 0; y < 800; y++) {
        if (tbl[y].xleft < tbl[y].xright) {
            DrawDDALine(hdc, {tbl[y].xleft, y}, {tbl[y].xright, y}, color);
        }
    }
}

void ConvexFill(HDC hdc, POINT p[], int n, COLORREF color) {
    EdgeTable tbl;
    init(tbl);
    polygon2table(p, n, tbl);
    table2screen(hdc, tbl, color);
}

void edge2tableGeneral(POINT &p1, POINT &p2, list<Node> tbl[]) {
    if (p1.y == p2.y) return;
    if (p1.y > p2.y) swap(p1, p2);
    tbl[p1.y].push_back(Node(p1.x, p2.y, (double)(p2.x - p1.x) / (p2.y - p1.y)));
}

void polygon2tableGeneral(POINT p[], int n, list<Node> tbl[]) {
    POINT point = p[n - 1];
    for (int i = 0; i < n; i++) {
        POINT nextPoint = p[i];
        edge2tableGeneral(point, nextPoint, tbl);

        point = p[i];
    }
}

void table2screenGeneral(HDC hdc, list<Node> tbl[], COLORREF color) {
    int y = 0;
    while (y < 80000 && tbl[y].empty()) y++;
    list<Node> active = tbl[y];
    while (!active.empty()) {
        active.sort();
        for (auto it = active.begin(); it != active.end(); ++it) {
            auto it2 = it++;
            int xleft = (int)ceil(it->x);
            int xright = (int)floor(it2->x);
            DrawDDALine(hdc, {xleft, y}, {xright, y}, color);
        }
        y++;
        for (auto it = active.begin(); it != active.end();) {
            if (it->ymax == y) {
                it = active.erase(it);
            } else {
                it->x += it->minv;
                ++it;
            }
        }
        for (auto &node : tbl[y]) {
            active.push_back(node);
        }
    }
}

void NonConvexFill(HDC hdc, POINT p[], int n, COLORREF color) {
    list<Node> tbl[80000];
    polygon2tableGeneral(p, n, tbl);
    table2screenGeneral(hdc, tbl, color);
}

POINT hermiteInterpolate(POINT p0, POINT m0, POINT p1, POINT m1, float t) {
    float t2 = t * t;
    float t3 = t2 * t;

    float h1 = 2 * t3 - 3 * t2 + 1;
    float h2 = t3 - 2 * t2 + t;
    float h3 = -2 * t3 + 3 * t2;
    float h4 = t3 - t2;

    POINT result;
    result.x = h1 * p0.x + h2 * m0.x + h3 * p1.x + h4 * m1.x;
    result.y = h1 * p0.y + h2 * m0.y + h3 * p1.y + h4 * m1.y;
    return result;
}

void drawCardinalSpline(HDC hdc, const vector<POINT>& points, float c, COLORREF color) {
    if (points.size() < 4) return;

    vector<POINT> q;
    q.push_back({0, 0});
    for (size_t i = 1; i < points.size() - 1; ++i) {
        int newx = (points[i + 1].x - points[i - 1].x) * (c / 2.0f);
        int newy = (points[i + 1].x - points[i - 1].x) * (c / 2.0f);
        q.push_back({newx, newy});
    }

    for (size_t i = 1; i < points.size() - 2; ++i) {
        for (float t = 0; t < 1.0f; t += 0.0001f) {
            POINT p0 = hermiteInterpolate(points[i], q[i], points[i + 1], q[i + 1], t);

            SetPixel(hdc, (int)p0.x, (int)p0.y, color);
        }
    }
}

void FillQuarterWithLines(HDC hdc, int xc, int yc, int R, int quarter, COLORREF c)
{
    for (int y = -R; y <= R; y++)
    {
        bool inSegment = false;
        POINT pStart, pEnd;

        for (int x = -R; x <= R; x++)
        {
            if (x * x + y * y <= R * R)
            {
                bool insideQuarter = false;

                switch (quarter)
                {
                    case 1: insideQuarter = (x >= 0 && y <= 0); break;
                    case 2: insideQuarter = (x <= 0 && y <= 0); break;
                    case 3: insideQuarter = (x <= 0 && y >= 0); break;
                    case 4: insideQuarter = (x >= 0 && y >= 0); break;
                }

                if (insideQuarter)
                {
                    if (!inSegment)
                    {
                        pStart = { xc + x, yc + y };
                        inSegment = true;
                    }
                    pEnd = { xc + x, yc + y };
                }
                else if (inSegment)
                {
                    DrawMidpointLine(hdc, pStart, pEnd, c);
                    inSegment = false;
                }
            }
            else if (inSegment)
            {
                DrawMidpointLine(hdc, pStart, pEnd, c);
                inSegment = false;
            }
        }

        if (inSegment)
        {
            DrawMidpointLine(hdc, pStart, pEnd, c);
        }
    }
}

void FillQuarterWithCircles(HDC hdc, int xc, int yc, int R, int quarter, COLORREF c)
{
    int smallR = 5;
    for (int y = -R; y <= R; y += 2 * smallR) {
        for (int x = -R; x <= R; x += 2 * smallR) {
            if (x * x + y * y <= R * R) {
                bool draw = false;
                if (quarter == 1 && x >= 0 && y <= 0) draw = true;
                else if (quarter == 2 && x <= 0 && y <= 0) draw = true;
                else if (quarter == 3 && x <= 0 && y >= 0) draw = true;
                else if (quarter == 4 && x >= 0 && y >= 0) draw = true;
                POINT p;
                p.x = xc + x;
                p.y = yc + y;
                if (draw) DrawModifiedMidpointCircle(hdc, p, smallR, c);
            }
        }
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
                pointCount = 0;
                currentAlgorithm = LOWORD(wParam);
                InitializeSquareClipWindow(hwnd);
                if (currentAlgorithm == ID_CLIP_SQUARE_POINT || currentAlgorithm == ID_CLIP_SQUARE_LINE) {
                    InitializeSquareClipWindow(hwnd);
                } else if (currentAlgorithm == ID_CLIP_RECT_POINT || currentAlgorithm == ID_CLIP_RECT_LINE) {
                    InitializeRectClipWindow(hwnd);
                }
                InvalidateRect(hwnd, NULL, TRUE); // Redraw to show clip window
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
        points[pointCount++] = pt;

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
            currentAlgorithm == ID_MOD_MIDPOINT_CIRCLE ||
            currentAlgorithm == ID_FILL_CIRCLE_CIRCLES ||
            currentAlgorithm == ID_FILL_CIRCLE_LINES)) {
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
                    case ID_FILL_CIRCLE_CIRCLES: {
                        DrawModifiedMidpointCircle(hdc, points[0], radius, currentColor);
                        int quarter;
                        cout << "Enter the quarter you wanna fill: ";
                        cin >> quarter;
                        if(quarter > 4 || quarter < 1) FillQuarterWithCircles(hdc, points[0].x, points[0].y, radius, 1, currentColor);
                        else FillQuarterWithCircles(hdc, points[0].x, points[0].y, radius, quarter, currentColor);
                        break;
                    }
                    case ID_FILL_CIRCLE_LINES: {
                        DrawModifiedMidpointCircle(hdc, points[0], radius, currentColor);
                        int quarter;
                        cout << "Enter the quarter you wanna fill: ";
                        cin >> quarter;
                        if(quarter > 4 || quarter < 1) FillQuarterWithLines(hdc, points[0].x, points[0].y, radius, 1, currentColor);
                        else FillQuarterWithLines(hdc, points[0].x, points[0].y, radius, quarter, currentColor);
                        break;
                    }
                }
                ReleaseDC(hwnd, hdc);
                pointCount = 0;
        } else if (pointCount == 3 &&
            (currentAlgorithm == ID_ELLIPSE_DIRECT ||
             currentAlgorithm == ID_ELLIPSE_POLAR ||
             currentAlgorithm == ID_ELLIPSE_MIDPOINT))
        {
            HDC hdc = GetDC(hwnd);
            int rx = abs(points[1].x - points[0].x);
            int ry = abs(points[2].y - points[0].y);
            switch (currentAlgorithm) {
                case ID_ELLIPSE_DIRECT:
                    DrawDirectEllipse(hdc, points[0], rx, ry, currentColor);
                    break;
                case ID_ELLIPSE_POLAR:
                    DrawPolarEllipse(hdc, points[0], rx, ry, currentColor);
                    break;
                case ID_ELLIPSE_MIDPOINT:
                    DrawMidpointEllipse(hdc, points[0], rx, ry, currentColor);
                    break;
            }
            ReleaseDC(hwnd, hdc);
            pointCount = 0;
        } else if (pointCount == 1 && (currentAlgorithm == ID_FLOOD_RECURSIVE || currentAlgorithm == ID_FLOOD_NONRECURSIVE)) {
            HDC hdc = GetDC(hwnd);
            switch (currentAlgorithm) {
                case ID_FLOOD_RECURSIVE:
                    FloodFillRecursive(hdc, points[0].x, points[0].y, bgColor, currentColor);
                    break;
                case ID_FLOOD_NONRECURSIVE:
                    FloodFillNonRecursive(hdc, points[0].x, points[0].y, bgColor, currentColor);
                    break;
            }
        } else if (pointCount == 2 && (currentAlgorithm == ID_FILL_SQUARE_HERMIT || currentAlgorithm == ID_FILL_RECTANGLE_BEZIER)) {
            HDC hdc = GetDC(hwnd);
            switch(currentAlgorithm) {
                case ID_FILL_SQUARE_HERMIT:
                    FillSquareWithVerticalHermite(hdc, points[0], points[1], currentColor);
                    break;
                case ID_FILL_RECTANGLE_BEZIER:
                    FillRectangleWithHorizontalBezier(hdc, points[0], points[1], currentColor);
                    break;
            }
            ReleaseDC(hwnd, hdc);
            pointCount = 0;
        } else if (pointCount == 1 && currentAlgorithm == ID_CLIP_SQUARE_POINT) {
            HDC hdc = GetDC(hwnd);
            bool clipped = ClipPointSquare(hdc, points[0], currentColor);
            if (clipped) {
                cout << "Point (" << points[0].x << ", " << points[0].y << ") is inside the clipping window" << endl;
            } else {
                cout << "Point (" << points[0].x << ", " << points[0].y << ") is clipped (outside window)" << endl;
            }
            ReleaseDC(hwnd, hdc);
            pointCount = 0;
        }
        else if (pointCount == 2 && currentAlgorithm == ID_CLIP_SQUARE_LINE) {
            HDC hdc = GetDC(hwnd);
            ClipLineSquare(hdc, points[0], points[1], currentColor);
            cout << "Line from (" << points[0].x << ", " << points[0].y << ") to ("
                 << points[1].x << ", " << points[1].y << ") processed for clipping" << endl;
            ReleaseDC(hwnd, hdc);
            pointCount = 0;
        } else if (pointCount == 1 && currentAlgorithm == ID_CLIP_RECT_POINT) {
            HDC hdc = GetDC(hwnd);
            bool clipped = ClipPointRect(hdc, points[0], currentColor);
            if (clipped) {
                cout << "Point (" << points[0].x << ", " << points[0].y << ") is inside the rectangle clipping window" << endl;
            } else {
                cout << "Point (" << points[0].x << ", " << points[0].y << ") is clipped (outside rectangle window)" << endl;
            }
            ReleaseDC(hwnd, hdc);
            pointCount = 0;
        }
        else if (pointCount == 2 && currentAlgorithm == ID_CLIP_RECT_LINE) {
            HDC hdc = GetDC(hwnd);
            ClipLineRect(hdc, points[0], points[1], currentColor);
            cout << "Line from (" << points[0].x << ", " << points[0].y << ") to ("
                 << points[1].x << ", " << points[1].y << ") processed for rectangle clipping" << endl;
            ReleaseDC(hwnd, hdc);
            pointCount = 0;
        }
        //TODO: Add the other algorithms.

        break;
    }
    case WM_RBUTTONDOWN: {
        if (pointCount > 2 && (currentAlgorithm == ID_FILL_CONVEX || currentAlgorithm == ID_FILL_NONCONVEX)) {
            HDC hdc = GetDC(hwnd);
            switch(currentAlgorithm) {
                case ID_FILL_CONVEX:
                    ConvexFill(hdc, points, pointCount, currentColor);
                    break;
                case ID_FILL_NONCONVEX:
                    NonConvexFill(hdc, points, pointCount, currentColor);
                    break;
            }
            ReleaseDC(hwnd, hdc);
            pointCount = 0;
        } else if (pointCount > 0 && currentAlgorithm == ID_CARDINAL_SPLINE) {
            HDC hdc = GetDC(hwnd);
            drawCardinalSpline(hdc, vector<POINT>(points, points + pointCount), 0.5f, currentColor);
            ReleaseDC(hwnd, hdc);
            pointCount = 0;
        }
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

        // Draw clipping window if using square clipping algorithms
        if (currentAlgorithm == ID_CLIP_SQUARE_POINT || currentAlgorithm == ID_CLIP_SQUARE_LINE) {
            DrawSquareClipWindow(hdc);
        } else if (currentAlgorithm == ID_CLIP_RECT_POINT || currentAlgorithm == ID_CLIP_RECT_LINE) {
            DrawRectClipWindow(hdc);
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
    AllocConsole();
    freopen("CONOUT$", "w", stdout);
    freopen("CONIN$", "r", stdin);
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
