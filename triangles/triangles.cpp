#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

//Tutaj wazne zmienne rozmiary ekranu.
const int canvas_width = 20;
const int canvas_height = 10;

class Line;
class Canvas;
class Point {
    public:
        Point(): x(0), y(0) {}
        Point(int _x, int _y): x(_x), y(_y) {}
        int GetX() const {return x;}
        int GetY() const {return y;}
        bool operator==(const Point &other) const {
            return x == other.GetX() && y == other.GetY();
        }
        friend Point operator+(const Point& p, const Line& vec);
    private:
        int x;
        int y;
};

class Line {
    public:
        Line(const Point& a, const Point& b): start(a), end(b) {}
        Line(): start(Point()), end(Point()){}
        void set_start(Point p) {start = p;}
        void set_stop(Point p) {end = p;}
        bool isPoint() const {return start == end;}
        Point GetStart() const {return start;}
        Point GetEnd() const {return end;}
    private:
        Point start;
        Point end;
};

Point operator+(const Point& p, const Line& vec) {
    auto t_x = vec.GetEnd().GetX() - vec.GetStart().GetX();
    auto t_y = vec.GetEnd().GetY() - vec.GetStart().GetY();
    return Point(p.x + t_x, p.y + t_y);
}

class Figure {
    public:
        Figure(): height(0), width(0), color('#'), base(Line()), pos(Point()) {}

        virtual float get_area() const = 0;
        virtual void draw(Canvas& c) = 0;

        virtual void move(Line vector) {
            base.set_start(base.GetStart() + vector);
            base.set_stop(base.GetEnd() + vector);
            pos = pos + vector;
        }

        virtual void set_height(int h) {
            height = h;
            RecalculateBase();
        }

        virtual void set_width(int w) {
            width = w;
            RecalculateBase();
        }

        unsigned int get_height() const {
            return height;
        }

        unsigned int get_width() const {
            return width;
        }

        void set_color(char c) {color = c;}

        void RecalculateBase() {
            Point basePositionA = Point(pos.GetX(), pos.GetY() + height);
            Point basePositionB = Point(pos.GetX()+width, pos.GetY() + height);
            base = Line(basePositionA, basePositionB);
        }

        Point get_position() const {
            return pos;
        }
    protected:
        unsigned int height;
        unsigned int width;
        char color;
        Line base;
        Point pos;
        
};

class Canvas {
    public:
        Canvas(): height(canvas_height), width(canvas_width), content_table(height, std::vector<char>(width, ' ')){}
        Canvas(int _height, int _width): height(_height), width(_width), content_table(height, std::vector<char>(width, ' ')) {}
        
        //Ja tego nie napisalem, mialem nie wydajna wersje ktora poprawil chat.
        std::string show() const {
            const size_t line_size = width + 1;  // znaki + \n
            const size_t total_size = height * line_size;
            std::string result(total_size, '\n');  // inicjalizuj wszystkimi \n
            char* ptr = result.data();
            for(const auto& line : content_table) {
                std::copy_n(line.data(), width, ptr);
                ptr += line_size;  // przeskakujemy do następnej linii
            }
            return result;
        }

        void Tick(int i, int j, char c) {
            if((j>=0 && j<height) && i>=0 && i<width) content_table[j][i] = c;
        }

        void DrawLine(Line l, char c) {
                auto dx = l.GetEnd().GetX() - l.GetStart().GetX();
                auto dy = l.GetEnd().GetY() - l.GetStart().GetY();
                std::abs(dx) == 0 ? dx = 0: dx = dx/std::abs(dx);
                std::abs(dy) == 0 ? dy = 0 : dy = dy/std::abs(dy);
                auto start_x = l.GetStart().GetX();
                auto start_y = l.GetStart().GetY();
                auto end = l.GetEnd();
                while(true) {
                    Tick(start_x, start_y, c);
                    if(start_x == end.GetX() && start_y == end.GetY()) break;
                    if(start_x != end.GetX()) start_x += dx;
                    if(start_y != end.GetY()) start_y += dy;
                }
        }

        private:    
            int height;
            int width;
            std::vector<std::vector<char>> content_table;
};
class Rectangle: public Figure {
    public:
        Rectangle(): Figure(){}
        void draw(Canvas& c) override {
            auto ceiling = GetCeiling();
            c.DrawLine(base,color);
            c.DrawLine(ceiling,color);
            c.DrawLine(Line(base.GetStart(),ceiling.GetStart()),color);
            c.DrawLine(Line(base.GetEnd(), ceiling.GetEnd()), color);
        }

        float get_area() const override {
            return width*height;
        }

        Line GetCeiling() {
            Point baseStart = base.GetStart();
            Point baseEnd = base.GetEnd();
            Point ceilingStart(baseStart.GetX(), baseStart.GetY() - height);
            Point ceilingEnd(baseEnd.GetX(), baseEnd.GetY() - height);
            return Line(ceilingStart,ceilingEnd);
        }
};
class Square: public Rectangle {
    public:
        Square(): Rectangle() {}
        void set_width(int w) override {
            width = w;
            height = w;
            RecalculateBase();
            
        }
        void set_height(int h) override {
            height = h; 
            width = h;
            RecalculateBase();
        }
};

class Triangle: public Figure {
    public:
        Triangle() : Figure(){}
        float get_area() const override {
            return (1.0/2)*width*height;
        }
};

class HalfSquareTriangle: public Triangle {
    public:
        HalfSquareTriangle(): Triangle() {}

        void draw(Canvas& c) override {
            c.DrawLine(base, color);
            Line leftSide(base.GetStart(), pos);
            Line rightSide(pos,base.GetEnd());
            c.DrawLine(leftSide, color);
            c.DrawLine(rightSide, color);
        }

        float get_area() const override {
            return (1.0/2)*height*width;
        }

        void set_width(int w) override {
            width = w;
            height = w;
            RecalculateBase();
        }

        void set_height(int h) override {
            height = h; 
            width = h;
            RecalculateBase();
        }
};

class QuarterSquareTriangle: public Triangle {
    public:
        QuarterSquareTriangle(): Triangle() {}
        void draw(Canvas& c) override {
            c.DrawLine(base,color);
            auto start = base.GetStart();
            auto end = base.GetEnd();
            int inc = 1;
            while(start.GetX() + inc <= end.GetX() - inc) {
                c.Tick(start.GetX() + inc, start.GetY()+inc, color);
                c.Tick(end.GetX() - inc, start.GetY() + inc, color);
                ++inc;
            }
        }

        void set_width(int w) override {
            width = w;
            height = w/2;
            RecalculateBase();
        }

        void set_height(int h) override {
            height = h; 
            width = 2*h;
            RecalculateBase();
        }

        void RecalculateBase() {
            Point basePositionA = Point(pos.GetX(), pos.GetY() );
            Point basePositionB = Point(pos.GetX()+width, pos.GetY());
            base = Line(basePositionA, basePositionB);
        }
};


int main() {
    return 0;
}