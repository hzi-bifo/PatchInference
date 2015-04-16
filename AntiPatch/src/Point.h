/* 
 * File:   Point.h
 * Author: l.Klesper
 *
 * Created on July 4, 2012, 1:08 PM
 */

#ifndef POINT_H
#define	POINT_H

/** class point
 *
 *	store a 3d point
 */
class Point {
private:
    float co_x, co_y, co_z;
public:

    Point() {};

    Point(float x, float y, float z) : co_x(x), co_y(y), co_z(z) {};

    Point(const Point &orig) : co_x(orig.co_x), co_y(orig.co_y), co_z(orig.co_z) {};

    virtual ~Point() {};

    void setPoint(float x, float y, float z) {
        co_x = x;
        co_y = y;
        co_z = z;
    }

    float getX() {
        return co_x;
    };

    float getY() {
        return co_y;
    };

    float getZ() {
        return co_z;
    };
};

#endif	/* POINT_H */

