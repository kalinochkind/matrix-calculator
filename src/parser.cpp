#include "parser.h"

int priority[128];
bool rightassoc[128];

const std::vector<std::pair<token_type, std::string> > splitExpression(const std::string &expr)
{
    priority[int('+')] = 0;
    priority[int('-')] = 0;
    priority[int('*')] = 1;
    priority[int('/')] = 1;
    priority[int('%')] = 1;
    priority[int('^')] = 2;
    priority[int('_')] = 3;
    priority[int('=')] = -1;
    rightassoc[int('^')] = true;
    rightassoc[int('_')] = true;
    std::string func, num, poly;
    std::vector<std::pair<token_type, std::string> > ans;
    token_type last = TOKEN_LEFTPAR;
    for(char i : expr + ' ')
    {
        if('A' <= i && i <= 'Z')
        {
            if(func.size())
                ans.push_back({TOKEN_FUNC, func});
            if(num.size())
                ans.push_back({TOKEN_NUMBER, num});
            num = func = "";
            if(last == TOKEN_MATRIX || last == TOKEN_NUMBER || last == TOKEN_RIGHTPAR)
                ans.push_back({TOKEN_OP, "*"});
            ans.push_back({TOKEN_MATRIX, std::string() + i});
            last = TOKEN_MATRIX;
        }
        else if(i == '+' || i == '*' || i == '^' || i == '(' || i == ')' || i == '=' || i == '/' || i == '%')
        {
            if(i == '/' && poly.size())
            {
                poly.push_back('/');
                continue;
            }
            if(func.size())
                ans.push_back({TOKEN_FUNC, func});
            if(num.size())
                ans.push_back({TOKEN_NUMBER, num});
            num = func = "";
            if(i == '(')
            {
                if(last == TOKEN_MATRIX || last == TOKEN_NUMBER || last == TOKEN_RIGHTPAR)
                    ans.push_back({TOKEN_OP, "*"});
                ans.push_back({TOKEN_LEFTPAR, ""});
            }
            else if(i == ')')
                ans.push_back({TOKEN_RIGHTPAR, ""});
            else
                ans.push_back({TOKEN_OP, std::string() + i});
            last = ans.back().first;
        }
        else if(i == '-')
        {
            if(poly.size())
            {
                poly.push_back('-');
                continue;
            }
            if(func.size())
                ans.push_back({TOKEN_FUNC, func});
            if(num.size())
                ans.push_back({TOKEN_NUMBER, num});
            num = func = "";
            if(last == TOKEN_LEFTPAR || last == TOKEN_OP || last == TOKEN_COMMA)
            {
                ans.push_back({TOKEN_OP, "_"});  // unary
            }
            else
            {
                ans.push_back({TOKEN_OP, "-"});
            }
            last = TOKEN_OP;
        }
        else if(('0' <= i && i <= '9') || i == '.')
        {
            if(poly.size())
            {
                poly.push_back(i);
                continue;
            }
            if(func.size())
                ans.push_back({TOKEN_FUNC, func});
            func = "";
            num.push_back(i);
            last = TOKEN_NUMBER;
        }
        else if('a' <= i && i <= 'z')
        {
            if(num.size())
                ans.push_back({TOKEN_NUMBER, num});
            num = "";
            func.push_back(i);
            last = TOKEN_FUNC;
        }
        else if(i == ',')
        {
            if(func.size())
                ans.push_back({TOKEN_FUNC, func});
            if(num.size())
                ans.push_back({TOKEN_NUMBER, num});
            num = func = "";
            ans.push_back({TOKEN_COMMA, ""});
            last = TOKEN_COMMA;
        }
        else if(i == '$')
        {
            if(func.size())
                ans.push_back({TOKEN_FUNC, func});
            if(num.size())
                ans.push_back({TOKEN_NUMBER, num});
            num = func = "";
            ans.push_back({TOKEN_DOLLAR, ""});
            last = TOKEN_DOLLAR;
        }
        else if(i == '"' || i == '\'')
        {
            if(poly.size())
            {
                if(poly[0] != i)
                {
                    ans.clear();
                    return ans;
                }
                ans.push_back({TOKEN_POLY, poly.substr(1)});
                poly = "";
                last = TOKEN_POLY;
            }
            else
            {
                poly = i;
                num = func = "";
            }
        }
        else
        {
            if(i == ' ' && poly.size())
            {
                poly.push_back(i);
                continue;
            }
            if(func.size())
                ans.push_back({TOKEN_FUNC, func});
            if(num.size())
                ans.push_back({TOKEN_NUMBER, num});
            num = func = "";
        }
    }
    return ans;
}
