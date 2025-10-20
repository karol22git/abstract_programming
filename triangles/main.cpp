#include <optional>
#include <vector>
#include <string>
const int canvas_width = 30;
const int canvas_height = 15;
class Canvas {
    public:
        Canvas(): height(canvas_height), width(canvas_width){}
        Canvas(int _height, int _width): height(_height), width(_width), content_table(height, std::vector<char>(width, ' ')) {}
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
        private:    
            int height;
            int width;
            std::vector<std::vector<char>> content_table;
};

class Line;
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
class Canvas;
Point operator+(const Point& p, const Line& vec) {
    auto t_x = vec.GetEnd().GetX() - vec.GetStart().GetX();
    auto t_y = vec.GetEnd().GetY() - vec.GetStart().GetY();
    return Point(p.x + t_x, p.y + t_y);
}
class Figure {
    public:
        virtual void move(Line vector) {
            base.set_start(base.GetStart() + vector);
            base.set_stop(base.GetStart() + vector);
            if(extra_point) extra_point.emplace(extra_point.value() + vector);
        }
        void set_color(char c) {color = c;}
        virtual void draw(Canvas* c) = 0;
    protected:
        unsigned int height;
        unsigned int width;
        char color = '#';
        Line base;
        std::optional<Point> extra_point = std::nullopt;
        
};

class Rectangle: public Figure {
    public:
        Rectangle(){}
        void draw(Canvas* c) override {

        }
};

int main() {

    return 0;
}