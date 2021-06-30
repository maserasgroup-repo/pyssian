
from pygments.style import Style
from pygments.token import (Keyword, Name, Comment, String, Error, 
                            Number, Operator, Generic, Literal,Text)

class OneDarkPro(Style):
    """
    Style adapted from: 
    https://github.com/Binaryify/OneDark-Pro
    """
    #default_style = 'OneDark-Pro'
    #title = ''

    background_color = "#282c34" 

    styles = {
        Comment.Multiline:              "italic #7f848e",
        Comment.Preproc:                "bold #7f848e",
        Comment.Single:                 "italic #7f848e",
        Comment.Special:                "bold italic #7f848e",
        Comment:                        "italic #7f848e",
        Error:                          "bg:#181a1f #c24038",
        Generic.Deleted:                "bg:#282c34 #e06c75",
        Generic.Emph:                   "italic #000000",
        Generic.Error:                  "#aa0000",
        Generic.Heading:                "#999999",
        Generic.Inserted:               "bg:#ddffdd #000000",
        Generic.Output:                 "#888888",
        Generic.Prompt:                 "#555555",
        Generic.Strong:                 "bold",
        Generic.Subheading:             "#aaaaaa",
        Generic.Traceback:              "#aa0000",
        Keyword.Constant:               "bold #e5c07b",
        Keyword.Declaration:            "bold #c678dd",
        Keyword.Namespace:              "bold #c678dd",
        Keyword.Pseudo:                 "bold #c678dd",
        Keyword.Reserved:               "bold #c678dd",
        Keyword.Type:                   "bold #445588",
        Keyword:                        "bold #c678dd",
        Literal.Number.Float:           "#d19a66",
        Literal.Number.Hex:             "#d19a66",
        Literal.Number.Integer.Long:    "#d19a66",
        Literal.Number.Integer:         "#d19a66",
        Literal.Number.Oct:             "#d19a66",
        Literal.Number:                 "#d19a66",
        Literal.String.Backtick:        "#98c379",
        Literal.String.Char:            "#98c379",
        Literal.String.Doc:             "#98c379",
        Literal.String.Double:          "#98c379",
        Literal.String.Escape:          "#56b6c2",
        Literal.String.Heredoc:         "#c678dd",
        Literal.String.Interpol:        "#e5c07b",
        Literal.String.Other:           "#61afef",
        Literal.String.Regex:           "#d19a66",
        Literal.String.Single:          "#98c379",
        Literal.String.Symbol:          "#98c379",
        Literal.String:                 "#c678dd",
        Name.Attribute:                 "#008080",
        Name.Builtin.Pseudo:            "#999999",
        Name.Builtin:                   "#50a6e7",
        Name.Class:                     "bold #e5c07b",
        Name.Constant:                  "#008080",
        Name.Decorator:                 "bold #61afef",
        Name.Entity:                    "#800080",
        Name.Exception:                 "bold #61afef",
        Name.Function:                  "bold #e5c07b",
        Name.Label:                     "bold #e5c07b",
        Name.Namespace:                 "#e5c07b",
        Name.Tag:                       "#000080",
        Name.Variable.Class:            "#e5c07b",
        Name.Variable.Global:           "#e5c07b",
        Name.Variable.Instance:         "#e06c75",
        Name.Variable:                  "#008080",
        Name:                           "#e06c75",
        Operator.Word:                  "bold #c678dd",
        Operator:                       "bold #56b6c2",
        Text.Whitespace:                "#bbbbbb",
        Text:                           "#fbfbfb",
    }
