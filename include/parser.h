#ifndef _PARSER_H
#define _PARSER_H

#include <vector>
#include <string>

enum token_type
{
    TOKEN_MATRIX, TOKEN_OP, TOKEN_FUNC, TOKEN_NUMBER, TOKEN_LEFTPAR, TOKEN_RIGHTPAR, TOKEN_COMMA, TOKEN_DOLLAR
};

extern int priority[128];
extern bool rightassoc[128];

const std::vector<std::pair<token_type, std::string> > splitExpression(const std::string &expr);
#endif
