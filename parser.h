#ifndef _PARSER_H
#define _PARSER_H

#include <vector>
#include "matrix.h"

enum token_type
{
    TOKEN_MATRIX, TOKEN_OP, TOKEN_FUNC, TOKEN_NUMBER, TOKEN_LEFTPAR, TOKEN_RIGHTPAR
};


const std::vector<std::pair<token_type, std::string> > splitExpression(const std::string &expr)
{
    std::string priority[128];
    priority['+'] = "0+";
    priority['*'] = "1*";
    priority['^'] = "2^";
    std::string func, num;
    std::vector<std::pair<token_type, std::string> > ans;
    for (char i : expr + ' ')
    {
        if ('A' <= i && i <= 'Z')
        {
            if (func.size())
                ans.push_back({TOKEN_FUNC, func});
            if (num.size())
                ans.push_back({TOKEN_NUMBER, num});
            num = func = "";
            ans.push_back({TOKEN_MATRIX, std::string() + i});

        }
        else if (i == '+' || i == '*' || i == '^' || i == '(' || i == ')')
        {
            if (func.size())
                ans.push_back({TOKEN_FUNC, func});
            if (num.size())
                ans.push_back({TOKEN_NUMBER, num});
            num = func = "";
            if (i == '(')
                ans.push_back({TOKEN_LEFTPAR, ""});
            else if (i == ')')
                ans.push_back({TOKEN_RIGHTPAR, ""});
            else
                ans.push_back({TOKEN_OP, priority[int(i)]});
        }
        else if (i == '-' || ('0' <= i && i <= '9'))
        {
            if (func.size())
                ans.push_back({TOKEN_FUNC, func});
            func = "";
            num.push_back(i);
        }
        else if ('a' <= i && i <= 'z')
        {
            if (num.size())
                ans.push_back({TOKEN_NUMBER, num});
            num = "";
            func.push_back(i);
        }
        else
        {
            if (func.size())
                ans.push_back({TOKEN_FUNC, func});
            if (num.size())
                ans.push_back({TOKEN_NUMBER, num});
            num = func = "";
        }
    }
    return ans;
}

#endif
