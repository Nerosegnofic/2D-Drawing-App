#include <windows.h>
#include <commdlg.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <list>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdio>

#define _USE_MATH_DEFINES

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

// Global variables
HBITMAP loadedBitmap = NULL;
COLORREF currentColor = RGB(255, 0, 0);
COLORREF bgColor = RGB(255, 255, 255);
POINT points[200];
int pointCount = 0;
int currentAlgorithm = ID_DDA_LINE;
HCURSOR currentCursor = LoadCursor(NULL, IDC_ARROW);

// Data structures
struct square_clip_window {
    int left, right, top, bottom;
    int size;
};

struct rect_clip_window {
    int left, right, top, bottom;
    int width, height;
};

struct Node {
    double x, minv;
    int ymax;
    Node(double _x = 0, int _ymax = 0, double _minv = 0) 
        : x(_x), minv(_minv), ymax(_ymax) {}

    bool operator<(const Node& other) const {
        return x < other.x;
    }
};

union Outcode {
    struct {
        unsigned left : 1;
        unsigned right : 1;
        unsigned top : 1;
        unsigned bottom : 1;
    };
    unsigned all : 4;
};

// Clipping windows
square_clip_window squareClipWindow = {0, 0, 0, 0, 100};
rect_clip_window rectClipWindow = {0, 0, 0, 0, 150, 100};

// Function prototypes
void DrawDDALine(HDC hdc, POINT p1, POINT p2, COLORREF color);
void DrawMidpointLine(HDC hdc, POINT p1, POINT p2, COLORREF color);
void DrawParametricLine(HDC hdc, POINT p1, POINT p2, COLORREF color);
void Draw8Points(HDC hdc, POINT center, int x, int y, COLORREF color);
void DrawDirectCircle(HDC hdc, POINT center, int radius, COLORREF color);
void DrawPolarCircle(HDC hdc, POINT center, int radius, COLORREF color);
void DrawIterativePolarCircle(HDC hdc, POINT center, int radius, COLORREF color);
void DrawMidpointCircle(HDC hdc, POINT center, int radius, COLORREF color);
void DrawModifiedMidpointCircle(HDC hdc, POINT center, int radius, COLORREF color);
void DrawDirectEllipse(HDC hdc, POINT center, int rx, int ry, COLORREF color);
void DrawPolarEllipse(HDC hdc, POINT center, int rx, int ry, COLORREF color);
void DrawMidpointEllipse(HDC hdc, POINT center, int rx, int ry, COLORREF color);
void FloodFillRecursive(HDC hdc, int x, int y, COLORREF targetColor, COLORREF fillColor);
void FloodFillNonRecursive(HDC hdc, int startX, int startY, COLORREF targetColor, COLORREF fillColor);
void FillSquareWithVerticalHermite(HDC hdc, POINT p1, POINT p2, COLORREF color);
void FillRectangleWithHorizontalBezier(HDC hdc, POINT p1, POINT p2, COLORREF color);
void ConvexFill(HDC hdc, POINT p[], int n, COLORREF color);
void NonConvexFill(HDC hdc, POINT p[], int n, COLORREF color);
void drawCardinalSpline(HDC hdc, const vector<POINT>& points, float c, COLORREF color);
void FillQuarterWithLines(HDC hdc, int xc, int yc, int R, int quarter, COLORREF c);
void FillQuarterWithCircles(HDC hdc, int xc, int yc, int R, int quarter, COLORREF c);
void InitializeSquareClipWindow(HWND hwnd);
void InitializeRectClipWindow(HWND hwnd);
void DrawSquareClipWindow(HDC hdc);
void DrawRectClipWindow(HDC hdc);
Outcode GetOutcode(POINT p, int xl, int xr, int yb, int yt);
POINT VIntersect(POINT p1, POINT p2, int x_edge);
POINT HIntersect(POINT p1, POINT p2, int y_edge);
bool ClipPointSquare(HDC hdc, POINT p, COLORREF color);
bool ClipPointRect(HDC hdc, POINT p, COLORREF color);
void ClipLineSquare(HDC hdc, POINT p1, POINT p2, COLORREF color);
void ClipLineRect(HDC hdc, POINT p1, POINT p2, COLORREF color);
void ClipPolygonRect(HDC hdc, POINT* vertices, int vertexCount, COLORREF color);
void SaveScreenToFile(HWND hwnd);
void LoadImageFromFile(HWND hwnd);
void AddMenus(HWND hwnd);
LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int);

// Line Drawing Algorithms
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
    } 
    else {
        if (p1.y > p2.y) swap(p1, p2);
        float mInv = (float)dx / dy;
        float x = p1.x;
        for (int y = p1.y; y <= p2.y; ++y) {
            SetPixel(hdc, round(x), y, color);
            x += mInv;
        }
    }
}

void DrawMidpointLine(HDC hdc, POINT p1, POINT p2, COLORREF color) {
    int dx = p2.x - p1.x;
    int dy = p2.y - p1.y;
    int x = p1.x, y = p1.y;

    int xInc = (dx > 0) ? 1 : -1;
    int yInc = (dy > 0) ? 1 : -1;
    dx = abs(dx);
    dy = abs(dy);

    if (dx >= dy) {
        int d = 2 * dy - dx;
        int d1 = 2 * dy;
        int d2 = 2 * (dy - dx);

        SetPixel(hdc, x, y, color);
        while (x != p2.x) {
            if (d < 0) d += d1;
            else {
                d += d2;
                y += yInc;
            }
            x += xInc;
            SetPixel(hdc, x, y, color);
        }
    } 
    else {
        int d = 2 * dx - dy;
        int d1 = 2 * dx;
        int d2 = 2 * (dx - dy);

        SetPixel(hdc, x, y, color);
        while (y != p2.y) {
            if (d < 0) d += d1;
            else {
                d += d2;
                x += xInc;
            }
            y += yInc;
            SetPixel(hdc, x, y, color);
        }
    }
}

void DrawParametricLine(HDC hdc, POINT p1, POINT p2, COLORREF color) {
    for (float t = 0; t <= 1; t += 0.001) {
        int x = p1.x + t * (p2.x - p1.x);
        int y = p1.y + t * (p2.y - p1.y);
        SetPixel(hdc, x, y, color);
    }
}

// Circle Drawing Algorithms
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

void DrawDirectCircle(HDC hdc, POINT center, int radius, COLORREF color) {
    int x = 0, y = radius;
    int radiusSq = radius * radius;

    while (x <= y) {
        Draw8Points(hdc, center, x, y, color);
        x++;
        y = round(sqrt(radiusSq - x * x));
    }
}

void DrawPolarCircle(HDC hdc, POINT center, int radius, COLORREF color) {
    double dtheta = 1.0 / radius;
    for (double theta = 0; theta <= M_PI/4; theta += dtheta) {
        int x = round(radius * cos(theta));
        int y = round(radius * sin(theta));
        Draw8Points(hdc, center, x, y, color);
    }
}

void DrawIterativePolarCircle(HDC hdc, POINT center, int radius, COLORREF color) {
    double x = radius, y = 0;
    double dtheta = 1.0 / radius;
    double cdtheta = cos(dtheta), sdtheta = sin(dtheta);

    while (x >= y) {
        Draw8Points(hdc, center, round(x), round(y), color);
        double x_new = x * cdtheta - y * sdtheta;
        y = x * sdtheta + y * cdtheta;
        x = x_new;
    }
}

void DrawMidpointCircle(HDC hdc, POINT center, int radius, COLORREF color) {
    int x = 0, y = radius;
    int d = 1 - radius;

    while (x <= y) {
        Draw8Points(hdc, center, x, y, color);
        if (d < 0) d += 2 * x + 3;
        else {
            d += 2 * (x - y) + 5;
            y--;
        }
        x++;
    }
}

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
        } 
        else {
            d += d2;
            d1 += 2;
            d2 += 4;
            y--;
        }
        x++;
    }
}

// Ellipse Drawing Algorithms
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

void DrawPolarEllipse(HDC hdc, POINT center, int rx, int ry, COLORREF color) {
    if (rx <= 0 || ry <= 0) return;
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

void DrawMidpointEllipse(HDC hdc, POINT center, int rx, int ry, COLORREF color) {
    if (rx <= 0 || ry <= 0) return;
    double rx2 = rx * rx;
    double ry2 = ry * ry;
    double x = 0, y = ry;
    double p = ry2 - rx2 * ry + 0.25 * rx2;

    while (ry2 * x < rx2 * y) {
        SetPixel(hdc, center.x + x, center.y + y, color);
        SetPixel(hdc, center.x - x, center.y + y, color);
        SetPixel(hdc, center.x + x, center.y - y, color);
        SetPixel(hdc, center.x - x, center.y - y, color);
        if (p < 0) p += 2 * ry2 * x + ry2;
        else {
            y--;
            p += 2 * ry2 * x - 2 * rx2 * y + ry2;
        }
        x++;
    }

    p = ry2 * (x + 0.5) * (x + 0.5) + rx2 * (y - 1) * (y - 1) - rx2 * ry2;
    while (y >= 0) {
        SetPixel(hdc, center.x + x, center.y + y, color);
        SetPixel(hdc, center.x - x, center.y + y, color);
        SetPixel(hdc, center.x + x, center.y - y, color);
        SetPixel(hdc, center.x - x, center.y - y, color);
        if (p > 0) p += -2 * rx2 * y + rx2;
        else {
            x++;
            p += 2 * ry2 * x - 2 * rx2 * y + rx2;
        }
        y--;
    }
}

// Filling Algorithms
void FloodFillRecursive(HDC hdc, int x, int y, COLORREF targetColor, COLORREF fillColor) {
    if (GetPixel(hdc, x, y) != targetColor || GetPixel(hdc, x, y) == fillColor) return;
    SetPixel(hdc, x, y, fillColor);
    FloodFillRecursive(hdc, x + 1, y, targetColor, fillColor);
    FloodFillRecursive(hdc, x - 1, y, targetColor, fillColor);
    FloodFillRecursive(hdc, x, y + 1, targetColor, fillColor);
    FloodFillRecursive(hdc, x, y - 1, targetColor, fillColor);
}

void FloodFillNonRecursive(HDC hdc, int startX, int startY, COLORREF targetColor, COLORREF fillColor) {
    stack<POINT> pixelStack;
    pixelStack.push({startX, startY});

    while (!pixelStack.empty()) {
        POINT current = pixelStack.top();
        pixelStack.pop();
        int x = current.x;
        int y = current.y;

        if (GetPixel(hdc, x, y) != targetColor || GetPixel(hdc, x, y) == fillColor) continue;
        SetPixel(hdc, x, y, fillColor);
        pixelStack.push({x + 1, y});
        pixelStack.push({x - 1, y});
        pixelStack.push({x, y + 1});
        pixelStack.push({x, y - 1});
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
    int totalSteps = 1000;
    vector<vector<int>> allY(endX - startX + 1);

    for (int xi = 0; xi <= endX - startX; xi++) {
        int x = startX + xi;
        vector<int>& pixelsY = allY[xi];
        double a0 = y1;
        double a1 = t0y;
        double a2 = -3 * y1 + 3 * y0 - 2 * t0y - t1y;
        double a3 = 2 * y1 - 2 * y0 + t0y + t1y;
        for (double t = 0; t <= 1.0; t += 0.001) {
            double y = a0 + a1 * t + a2 * t * t + a3 * t * t * t;
            pixelsY.push_back((int)y);
        }
    }

    for (int step = 0; step < totalSteps; step++) {
        for (int xi = 0; xi <= endX - startX; xi++) {
            if (step < allY[xi].size()) {
                SetPixel(hdc, startX + xi, allY[xi][step], color);
            }
        }
    }
}

void FillRectangleWithHorizontalBezier(HDC hdc, POINT p1, POINT p2, COLORREF color) {
    int width = abs(p2.x - p1.x);
    int height = abs(p2.y - p1.y);
    int dirX = (p2.x > p1.x) ? 1 : -1;
    int dirY = (p2.y > p1.y) ? 1 : -1;
    int xStart = p1.x;
    int xEnd = p1.x + dirX * width;

    HPEN pen = CreatePen(PS_SOLID, 1, color);
    HPEN oldPen = (HPEN)SelectObject(hdc, pen);

    for (int x = min(xStart, xEnd); x <= max(xStart, xEnd); x++) {
        POINT top = {x, p1.y};
        POINT ctrl1 = {x, p1.y + dirY * height / 3};
        POINT ctrl2 = {x, p1.y + dirY * 2 * height / 3};
        POINT bottom = {x, p1.y + dirY * height};
        
        MoveToEx(hdc, x, top.y, NULL);
        for (double t = 0.0; t <= 1.0; t += 0.01) {
            double y = pow(1 - t, 3) * top.y +
                       3 * pow(1 - t, 2) * t * ctrl1.y +
                       3 * (1 - t) * t * t * ctrl2.y +
                       t * t * t * bottom.y;
            LineTo(hdc, x, (int)y);
        }
    }
    SelectObject(hdc, oldPen);
    DeleteObject(pen);
}

typedef struct {
    int xleft, xright;
} EdgeTable[80000];

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
    double d = (double)(p2.x - p1.x) / (p2.y - p1.y);
    while (y < p2.y) {
        if (x < tbl[y].xleft) tbl[y].xleft = (int)ceil(x);
        if (x > tbl[y].xright) tbl[y].xright = (int)floor(x);
        y++;
        x += d;
    }
}

void polygon2table(POINT p[], int n, EdgeTable tbl) {
    POINT point = p[n - 1];
    for (int i = 0; i < n; i++) {
        edge2table(point, p[i], tbl);
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
        edge2tableGeneral(point, p[i], tbl);
        point = p[i];
    }
}

void table2screenGeneral(HDC hdc, list<Node> tbl[], COLORREF color) {
    int y = 0;
    while (y < 80000 && tbl[y].empty()) y++;
    list<Node> active = tbl[y];
    
    while (!active.empty()) {
        active.sort();
        auto it = active.begin();
        while (it != active.end()) {
            auto next = it;
            advance(next, 1);
            if (next == active.end()) break;
            
            int xleft = (int)ceil(it->x);
            int xright = (int)floor(next->x);
            if (xleft < xright) {
                DrawDDALine(hdc, {xleft, y}, {xright, y}, color);
            }
            advance(it, 2);
        }
        
        y++;
        for (auto it = active.begin(); it != active.end(); ) {
            if (it->ymax == y) it = active.erase(it);
            else {
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

void FillQuarterWithLines(HDC hdc, int xc, int yc, int R, int quarter, COLORREF c) {
    for (int y = -R; y <= R; y++) {
        bool inSegment = false;
        POINT pStart, pEnd;

        for (int x = -R; x <= R; x++) {
            if (x * x + y * y <= R * R) {
                bool insideQuarter = false;
                switch (quarter) {
                    case 1: insideQuarter = (x >= 0 && y <= 0); break;
                    case 2: insideQuarter = (x <= 0 && y <= 0); break;
                    case 3: insideQuarter = (x <= 0 && y >= 0); break;
                    case 4: insideQuarter = (x >= 0 && y >= 0); break;
                }

                if (insideQuarter) {
                    if (!inSegment) {
                        pStart = {xc + x, yc + y};
                        inSegment = true;
                    }
                    pEnd = {xc + x, yc + y};
                }
                else if (inSegment) {
                    DrawMidpointLine(hdc, pStart, pEnd, c);
                    inSegment = false;
                }
            }
            else if (inSegment) {
                DrawMidpointLine(hdc, pStart, pEnd, c);
                inSegment = false;
            }
        }
        if (inSegment) DrawMidpointLine(hdc, pStart, pEnd, c);
    }
}

void FillQuarterWithCircles(HDC hdc, int xc, int yc, int R, int quarter, COLORREF c) {
    int smallR = 5;
    for (int y = -R; y <= R; y += 2 * smallR) {
        for (int x = -R; x <= R; x += 2 * smallR) {
            if (x * x + y * y <= R * R) {
                bool draw = false;
                if (quarter == 1 && x >= 0 && y <= 0) draw = true;
                else if (quarter == 2 && x <= 0 && y <= 0) draw = true;
                else if (quarter == 3 && x <= 0 && y >= 0) draw = true;
                else if (quarter == 4 && x >= 0 && y >= 0) draw = true;
                
                if (draw) {
                    POINT p = {xc + x, yc + y};
                    DrawModifiedMidpointCircle(hdc, p, smallR, c);
                }
            }
        }
    }
}

POINT hermiteInterpolate(POINT p0, POINT m0, POINT p1, POINT m1, float t) {
    float t2 = t * t;
    float t3 = t2 * t;
    float h1 = 2 * t3 - 3 * t2 + 1;
    float h2 = t3 - 2 * t2 + t;
    float h3 = -2 * t3 + 3 * t2;
    float h4 = t3 - t2;

    return {
        (LONG)(h1 * p0.x + h2 * m0.x + h3 * p1.x + h4 * m1.x),
        (LONG)(h1 * p0.y + h2 * m0.y + h3 * p1.y + h4 * m1.y)
    };
}

void drawCardinalSpline(HDC hdc, const vector<POINT>& points, float c, COLORREF color) {
    if (points.size() < 4) return;
    vector<POINT> tangents;
    tangents.push_back({0, 0});
    
    for (size_t i = 1; i < points.size() - 1; ++i) {
        tangents.push_back({
            (LONG)((points[i + 1].x - points[i - 1].x) * (c / 2.0f)),
            (LONG)((points[i + 1].y - points[i - 1].y) * (c / 2.0f))
        });
    }

    for (size_t i = 1; i < points.size() - 2; ++i) {
        for (float t = 0; t < 1.0f; t += 0.001f) {
            POINT p = hermiteInterpolate(points[i], tangents[i], points[i+1], tangents[i+1], t);
            SetPixel(hdc, p.x, p.y, color);
        }
    }
}

// Clipping Functions
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

void DrawSquareClipWindow(HDC hdc) {
    HPEN pen = CreatePen(PS_SOLID, 2, RGB(255, 0, 0));
    HPEN oldPen = (HPEN)SelectObject(hdc, pen);
    Rectangle(hdc, squareClipWindow.left, squareClipWindow.top, 
              squareClipWindow.right + 1, squareClipWindow.bottom + 1);
    SelectObject(hdc, oldPen);
    DeleteObject(pen);
}

void DrawRectClipWindow(HDC hdc) {
    HPEN pen = CreatePen(PS_SOLID, 2, RGB(0, 0, 255));
    HPEN oldPen = (HPEN)SelectObject(hdc, pen);
    Rectangle(hdc, rectClipWindow.left, rectClipWindow.top,
              rectClipWindow.right + 1, rectClipWindow.bottom + 1);
    SelectObject(hdc, oldPen);
    DeleteObject(pen);
}

Outcode GetOutcode(POINT p, int xl, int xr, int yb, int yt) {
    Outcode result = {0};
    if (p.x < xl) result.left = 1;
    if (p.x > xr) result.right = 1;
    if (p.y < yb) result.bottom = 1;
    if (p.y > yt) result.top = 1;
    return result;
}

POINT VIntersect(POINT p1, POINT p2, int x_edge) {
    return {
        x_edge,
        (LONG)(p1.y + (x_edge - p1.x) * (double)(p2.y - p1.y) / (p2.x - p1.x))
    };
}

POINT HIntersect(POINT p1, POINT p2, int y_edge) {
    return {
        (LONG)(p1.x + (y_edge - p1.y) * (double)(p2.x - p1.x) / (p2.y - p1.y)),
        y_edge
    };
}

bool ClipPointSquare(HDC hdc, POINT p, COLORREF color) {
    if (p.x >= squareClipWindow.left && p.x <= squareClipWindow.right &&
        p.y >= squareClipWindow.top && p.y <= squareClipWindow.bottom) {
        SetPixel(hdc, p.x, p.y, color);
        return true;
    }
    return false;
}

bool ClipPointRect(HDC hdc, POINT p, COLORREF color) {
    if (p.x >= rectClipWindow.left && p.x <= rectClipWindow.right &&
        p.y >= rectClipWindow.top && p.y <= rectClipWindow.bottom) {
        SetPixel(hdc, p.x, p.y, color);
        return true;
    }
    return false;
}

void ClipLineSquare(HDC hdc, POINT p1, POINT p2, COLORREF color) {
    Outcode out1 = GetOutcode(p1, squareClipWindow.left, squareClipWindow.right, 
                              squareClipWindow.top, squareClipWindow.bottom);
    Outcode out2 = GetOutcode(p2, squareClipWindow.left, squareClipWindow.right, 
                              squareClipWindow.top, squareClipWindow.bottom);

    while (true) {
        if (out1.all == 0 && out2.all == 0) {
            DrawDDALine(hdc, p1, p2, color);
            return;
        }
        if (out1.all & out2.all) return;

        Outcode outcode = out1.all != 0 ? out1 : out2;
        POINT* point = out1.all != 0 ? &p1 : &p2;

        if (outcode.left) *point = VIntersect(p1, p2, squareClipWindow.left);
        else if (outcode.right) *point = VIntersect(p1, p2, squareClipWindow.right);
        else if (outcode.bottom) *point = HIntersect(p1, p2, squareClipWindow.top);
        else if (outcode.top) *point = HIntersect(p1, p2, squareClipWindow.bottom);

        if (point == &p1) out1 = GetOutcode(p1, squareClipWindow.left, squareClipWindow.right, 
                                           squareClipWindow.top, squareClipWindow.bottom);
        else out2 = GetOutcode(p2, squareClipWindow.left, squareClipWindow.right, 
                              squareClipWindow.top, squareClipWindow.bottom);
    }
}

void ClipLineRect(HDC hdc, POINT p1, POINT p2, COLORREF color) {
    Outcode out1 = GetOutcode(p1, rectClipWindow.left, rectClipWindow.right,
                             rectClipWindow.top, rectClipWindow.bottom);
    Outcode out2 = GetOutcode(p2, rectClipWindow.left, rectClipWindow.right,
                             rectClipWindow.top, rectClipWindow.bottom);

    while (true) {
        if (out1.all == 0 && out2.all == 0) {
            DrawDDALine(hdc, p1, p2, color);
            return;
        }
        if (out1.all & out2.all) return;

        Outcode outcode = out1.all != 0 ? out1 : out2;
        POINT* point = out1.all != 0 ? &p1 : &p2;

        if (outcode.left) *point = VIntersect(p1, p2, rectClipWindow.left);
        else if (outcode.right) *point = VIntersect(p1, p2, rectClipWindow.right);
        else if (outcode.bottom) *point = HIntersect(p1, p2, rectClipWindow.top);
        else if (outcode.top) *point = HIntersect(p1, p2, rectClipWindow.bottom);

        if (point == &p1) out1 = GetOutcode(p1, rectClipWindow.left, rectClipWindow.right,
                                           rectClipWindow.top, rectClipWindow.bottom);
        else out2 = GetOutcode(p2, rectClipWindow.left, rectClipWindow.right,
                              rectClipWindow.top, rectClipWindow.bottom);
    }
}

vector<POINT> ClipPolygonAgainstEdge(const vector<POINT>& inputVertices, int edgeType, int edgeValue) {
    vector<POINT> outputVertices;
    if (inputVertices.empty()) return outputVertices;

    POINT prevVertex = inputVertices.back();
    for (const POINT& currentVertex : inputVertices) {
        bool currentInside = false;
        bool prevInside = false;

        switch (edgeType) {
            case 0: currentInside = (currentVertex.x >= edgeValue); 
                    prevInside = (prevVertex.x >= edgeValue); break;
            case 1: currentInside = (currentVertex.x <= edgeValue); 
                    prevInside = (prevVertex.x <= edgeValue); break;
            case 2: currentInside = (currentVertex.y >= edgeValue); 
                    prevInside = (prevVertex.y >= edgeValue); break;
            case 3: currentInside = (currentVertex.y <= edgeValue); 
                    prevInside = (prevVertex.y <= edgeValue); break;
        }

        if (currentInside) {
            if (!prevInside) {
                POINT intersection = (edgeType < 2) ? 
                    VIntersect(prevVertex, currentVertex, edgeValue) :
                    HIntersect(prevVertex, currentVertex, edgeValue);
                outputVertices.push_back(intersection);
            }
            outputVertices.push_back(currentVertex);
        } 
        else if (prevInside) {
            POINT intersection = (edgeType < 2) ? 
                VIntersect(prevVertex, currentVertex, edgeValue) :
                HIntersect(prevVertex, currentVertex, edgeValue);
            outputVertices.push_back(intersection);
        }
        prevVertex = currentVertex;
    }
    return outputVertices;
}

void ClipPolygonRect(HDC hdc, POINT* vertices, int vertexCount, COLORREF color) {
    if (vertexCount < 3) return;
    vector<POINT> clippedVertices(vertices, vertices + vertexCount);

    clippedVertices = ClipPolygonAgainstEdge(clippedVertices, 0, rectClipWindow.left);
    clippedVertices = ClipPolygonAgainstEdge(clippedVertices, 1, rectClipWindow.right);
    clippedVertices = ClipPolygonAgainstEdge(clippedVertices, 2, rectClipWindow.top);
    clippedVertices = ClipPolygonAgainstEdge(clippedVertices, 3, rectClipWindow.bottom);

    if (clippedVertices.size() >= 3) {
        for (size_t i = 0; i < clippedVertices.size(); i++) {
            POINT p1 = clippedVertices[i];
            POINT p2 = clippedVertices[(i + 1) % clippedVertices.size()];
            DrawDDALine(hdc, p1, p2, color);
        }
    }
}

// File Operations
void SaveScreenToFile(HWND hwnd) {
    HDC hdcScreen = GetDC(hwnd);
    HDC hdcMem = CreateCompatibleDC(hdcScreen);
    RECT rc;
    GetClientRect(hwnd, &rc);

    HBITMAP hBitmap = CreateCompatibleBitmap(hdcScreen, rc.right, rc.bottom);
    SelectObject(hdcMem, hBitmap);
    BitBlt(hdcMem, 0, 0, rc.right, rc.bottom, hdcScreen, 0, 0, SRCCOPY);

    char fileName[MAX_PATH] = "untitled.bmp";
    OPENFILENAME ofn = { sizeof(ofn) };
    ofn.lpstrFilter = "Bitmap Files (*.bmp)\0*.bmp\0All Files (*.*)\0*.*\0";
    ofn.lpstrFile = fileName;
    ofn.nMaxFile = MAX_PATH;
    ofn.lpstrDefExt = "bmp";
    ofn.lpstrTitle = "Save As Bitmap";
    ofn.Flags = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST;

    if (GetSaveFileName(&ofn)) {
        BITMAPFILEHEADER bmfHeader = {0};
        BITMAPINFOHEADER bi = {0};
        bi.biSize = sizeof(BITMAPINFOHEADER);
        bi.biWidth = rc.right;
        bi.biHeight = -rc.bottom;
        bi.biPlanes = 1;
        bi.biBitCount = 24;
        bi.biCompression = BI_RGB;

        DWORD dwBmpSize = ((rc.right * 3 + 3) & ~3) * abs(rc.bottom);
        BYTE* lpBits = new BYTE[dwBmpSize];
        GetDIBits(hdcMem, hBitmap, 0, rc.bottom, lpBits, (BITMAPINFO*)&bi, DIB_RGB_COLORS);

        ofstream file(fileName, ios::binary);
        if (file) {
            bmfHeader.bfType = 0x4D42;
            bmfHeader.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
            bmfHeader.bfSize = bmfHeader.bfOffBits + dwBmpSize;

            file.write((char*)&bmfHeader, sizeof(bmfHeader));
            file.write((char*)&bi, sizeof(bi));
            file.write((char*)lpBits, dwBmpSize);
        }
        delete[] lpBits;
    }

    DeleteObject(hBitmap);
    DeleteDC(hdcMem);
    ReleaseDC(hwnd, hdcScreen);
}

void LoadImageFromFile(HWND hwnd) {
    char fileName[MAX_PATH] = "";
    OPENFILENAME ofn = { sizeof(ofn) };
    ofn.hwndOwner = hwnd;
    ofn.lpstrFile = fileName;
    ofn.nMaxFile = MAX_PATH;
    ofn.lpstrFilter = "Bitmap Files (*.bmp)\0*.bmp\0All Files (*.*)\0*.*\0";
    ofn.lpstrTitle = "Select Bitmap File to Open";
    ofn.Flags = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST;

    if (GetOpenFileName(&ofn)) {
        HBITMAP hBitmap = (HBITMAP)LoadImage(
            NULL, fileName, IMAGE_BITMAP, 0, 0, 
            LR_LOADFROMFILE | LR_CREATEDIBSECTION
        );

        if (hBitmap) {
            if (loadedBitmap) DeleteObject(loadedBitmap);
            loadedBitmap = hBitmap;
            InvalidateRect(hwnd, NULL, TRUE);
        }
    }
}

// Menu System
void AddMenus(HWND hwnd) {
    HMENU hMain = CreateMenu();
    HMENU hColor = CreateMenu(), hCursor = CreateMenu(), hLines = CreateMenu(), hCircles = CreateMenu();
    HMENU hFills = CreateMenu(), hClipping = CreateMenu(), hEllipse = CreateMenu();

    AppendMenu(hColor, MF_STRING, ID_COLOR_RED, "Red");
    AppendMenu(hColor, MF_STRING, ID_COLOR_GREEN, "Green");
    AppendMenu(hColor, MF_STRING, ID_COLOR_BLUE, "Blue");

    AppendMenu(hCursor, MF_STRING, ID_CURSOR_ARROW, "Arrow");
    AppendMenu(hCursor, MF_STRING, ID_CURSOR_CROSS, "Cross");

    AppendMenu(hLines, MF_STRING, ID_DDA_LINE, "DDA Line");
    AppendMenu(hLines, MF_STRING, ID_MIDPOINT_LINE, "Midpoint Line");
    AppendMenu(hLines, MF_STRING, ID_PARAMETRIC_LINE, "Parametric Line");

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

// Window Procedure
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
                    
                case ID_CLEAR_SCREEN:
                    pointCount = 0;
                    if (loadedBitmap) {
                        DeleteObject(loadedBitmap);
                        loadedBitmap = NULL;
                    }
                    InvalidateRect(hwnd, NULL, TRUE);
                    break;
                    
                case ID_SAVE_IMAGE:
                    SaveScreenToFile(hwnd);
                    break;
                    
                case ID_LOAD_IMAGE:
                    LoadImageFromFile(hwnd);
                    break;
                    
                default: {
                    int previousAlgorithm = currentAlgorithm;
                    currentAlgorithm = LOWORD(wParam);
                    pointCount = 0;
                    
                    if (previousAlgorithm == ID_CLIP_SQUARE_POINT || 
                        previousAlgorithm == ID_CLIP_SQUARE_LINE ||
                        previousAlgorithm == ID_CLIP_RECT_POINT || 
                        previousAlgorithm == ID_CLIP_RECT_LINE || 
                        previousAlgorithm == ID_CLIP_RECT_POLYGON) {
                        InvalidateRect(hwnd, nullptr, TRUE);
                    }
                    
                    if (currentAlgorithm == ID_CLIP_SQUARE_POINT || 
                        currentAlgorithm == ID_CLIP_SQUARE_LINE) {
                        InitializeSquareClipWindow(hwnd);
                    } 
                    else if (currentAlgorithm == ID_CLIP_RECT_POINT || 
                            currentAlgorithm == ID_CLIP_RECT_LINE || 
                            currentAlgorithm == ID_CLIP_RECT_POLYGON) {
                        InitializeRectClipWindow(hwnd);
                    }
                    break;
                }
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
            HDC hdc = GetDC(hwnd);
            
            switch (currentAlgorithm) {
                case ID_DDA_LINE:
                case ID_MIDPOINT_LINE:
                case ID_PARAMETRIC_LINE:
                    if (pointCount == 2) {
                        if (currentAlgorithm == ID_DDA_LINE) DrawDDALine(hdc, points[0], points[1], currentColor);
                        else if (currentAlgorithm == ID_MIDPOINT_LINE) DrawMidpointLine(hdc, points[0], points[1], currentColor);
                        else DrawParametricLine(hdc, points[0], points[1], currentColor);
                        pointCount = 0;
                    }
                    break;
                    
                case ID_DIRECT_CIRCLE:
                case ID_POLAR_CIRCLE:
                case ID_ITER_POLAR_CIRCLE:
                case ID_MIDPOINT_CIRCLE:
                case ID_MOD_MIDPOINT_CIRCLE:
                case ID_FILL_CIRCLE_CIRCLES:
                case ID_FILL_CIRCLE_LINES:
                    if (pointCount == 2) {
                        int radius = static_cast<int>(sqrt(
                            pow(points[1].x - points[0].x, 2) +
                            pow(points[1].y - points[0].y, 2)
                        );
                        
                        if (currentAlgorithm == ID_DIRECT_CIRCLE) DrawDirectCircle(hdc, points[0], radius, currentColor);
                        else if (currentAlgorithm == ID_POLAR_CIRCLE) DrawPolarCircle(hdc, points[0], radius, currentColor);
                        else if (currentAlgorithm == ID_ITER_POLAR_CIRCLE) DrawIterativePolarCircle(hdc, points[0], radius, currentColor);
                        else if (currentAlgorithm == ID_MIDPOINT_CIRCLE) DrawMidpointCircle(hdc, points[0], radius, currentColor);
                        else if (currentAlgorithm == ID_MOD_MIDPOINT_CIRCLE) DrawModifiedMidpointCircle(hdc, points[0], radius, currentColor);
                        else if (currentAlgorithm == ID_FILL_CIRCLE_CIRCLES) {
                            DrawModifiedMidpointCircle(hdc, points[0], radius, currentColor);
                            int quarter;
                            cout << "Enter quarter (1-4): ";
                            cin >> quarter;
                            FillQuarterWithCircles(hdc, points[0].x, points[0].y, radius, quarter, currentColor);
                        }
                        else if (currentAlgorithm == ID_FILL_CIRCLE_LINES) {
                            DrawModifiedMidpointCircle(hdc, points[0], radius, currentColor);
                            int quarter;
                            cout << "Enter quarter (1-4): ";
                            cin >> quarter;
                            FillQuarterWithLines(hdc, points[0].x, points[0].y, radius, quarter, currentColor);
                        }
                        pointCount = 0;
                    }
                    break;
                    
                case ID_ELLIPSE_DIRECT:
                case ID_ELLIPSE_POLAR:
                case ID_ELLIPSE_MIDPOINT:
                    if (pointCount == 3) {
                        int rx = abs(points[1].x - points[0].x);
                        int ry = abs(points[2].y - points[0].y);
                        if (currentAlgorithm == ID_ELLIPSE_DIRECT) DrawDirectEllipse(hdc, points[0], rx, ry, currentColor);
                        else if (currentAlgorithm == ID_ELLIPSE_POLAR) DrawPolarEllipse(hdc, points[0], rx, ry, currentColor);
                        else DrawMidpointEllipse(hdc, points[0], rx, ry, currentColor);
                        pointCount = 0;
                    }
                    break;
                    
                case ID_FLOOD_RECURSIVE:
                case ID_FLOOD_NONRECURSIVE:
                    if (pointCount == 1) {
                        if (currentAlgorithm == ID_FLOOD_RECURSIVE) 
                            FloodFillRecursive(hdc, points[0].x, points[0].y, bgColor, currentColor);
                        else 
                            FloodFillNonRecursive(hdc, points[0].x, points[0].y, bgColor, currentColor);
                    }
                    break;
                    
                case ID_FILL_SQUARE_HERMIT:
                case ID_FILL_RECTANGLE_BEZIER:
                    if (pointCount == 2) {
                        if (currentAlgorithm == ID_FILL_SQUARE_HERMIT) 
                            FillSquareWithVerticalHermite(hdc, points[0], points[1], currentColor);
                        else 
                            FillRectangleWithHorizontalBezier(hdc, points[0], points[1], currentColor);
                        pointCount = 0;
                    }
                    break;
                    
                case ID_CLIP_SQUARE_POINT:
                    if (pointCount == 1) {
                        ClipPointSquare(hdc, points[0], currentColor);
                        pointCount = 0;
                    }
                    break;
                    
                case ID_CLIP_SQUARE_LINE:
                    if (pointCount == 2) {
                        ClipLineSquare(hdc, points[0], points[1], currentColor);
                        pointCount = 0;
                    }
                    break;
                    
                case ID_CLIP_RECT_POINT:
                    if (pointCount == 1) {
                        ClipPointRect(hdc, points[0], currentColor);
                        pointCount = 0;
                    }
                    break;
                    
                case ID_CLIP_RECT_LINE:
                    if (pointCount == 2) {
                        ClipLineRect(hdc, points[0], points[1], currentColor);
                        pointCount = 0;
                    }
                    break;
            }
            ReleaseDC(hwnd, hdc);
            break;
        }
            
        case WM_RBUTTONDOWN:
            if (pointCount > 0) {
                HDC hdc = GetDC(hwnd);
                switch (currentAlgorithm) {
                    case ID_FILL_CONVEX:
                    case ID_FILL_NONCONVEX:
                        if (pointCount > 2) {
                            if (currentAlgorithm == ID_FILL_CONVEX) 
                                ConvexFill(hdc, points, pointCount, currentColor);
                            else 
                                NonConvexFill(hdc, points, pointCount, currentColor);
                            pointCount = 0;
                        }
                        break;
                        
                    case ID_CARDINAL_SPLINE:
                        drawCardinalSpline(hdc, vector<POINT>(points, points + pointCount), 0.5f, currentColor);
                        pointCount = 0;
                        break;
                        
                    case ID_CLIP_RECT_POLYGON:
                        if (pointCount >= 3) {
                            HPEN grayPen = CreatePen(PS_SOLID, 1, RGB(192, 192, 192));
                            HPEN oldPen = (HPEN)SelectObject(hdc, grayPen);
                            for (int i = 0; i < pointCount; i++) {
                                POINT p1 = points[i];
                                POINT p2 = points[(i + 1) % pointCount];
                                DrawDDALine(hdc, p1, p2, RGB(192, 192, 192));
                            }
                            SelectObject(hdc, oldPen);
                            DeleteObject(grayPen);
                            ClipPolygonRect(hdc, points, pointCount, currentColor);
                            pointCount = 0;
                        }
                        break;
                }
                ReleaseDC(hwnd, hdc);
            }
            break;
            
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
            
            if (currentAlgorithm == ID_CLIP_SQUARE_POINT || currentAlgorithm == ID_CLIP_SQUARE_LINE) {
                DrawSquareClipWindow(hdc);
            } 
            else if (currentAlgorithm == ID_CLIP_RECT_POINT || 
                    currentAlgorithm == ID_CLIP_RECT_LINE || 
                    currentAlgorithm == ID_CLIP_RECT_POLYGON) {
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

// Entry Point
int WINAPI WinMain(HINSTANCE hInst, HINSTANCE, LPSTR, int nCmdShow) {
    AllocConsole();
    freopen("CONOUT$", "w", stdout);
    freopen("CONIN$", "r", stdin);
    
    WNDCLASS wc = {0};
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
