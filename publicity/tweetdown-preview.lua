-- Recover asterisks since Twitter italics are poor form
function Emph(elem)
  content_str = pandoc.utils.stringify(elem.c)
  return pandoc.Str("*" .. content_str .. "*")
end

-- Convert **strong** elements to ALL CAPS
function Strong(elem)
  content_str = pandoc.utils.stringify(elem.c)
  return pandoc.Str(string.upper(content_str))
end

wordcount = { 
  Str = function(el) 
    -- we don't count a word if it's entirely punctuation: 
    if el.text:match("%P") then 
        words = words + 1 
    end 
    -- URLs get shortened to 23 characters on Twitter
    if el.text:match("^http[s]*://") then
      charcount = charcount + 23
      print(el.text)
    else
      charcount = charcount + utf8.len(el.text)
    end
  end, 

  -- TODO: Which whitespace chars are (& are not) matched by Space?
  --       Does a LineBreak ever end up in the knitted md file?
  Space = function(el)
    charcount = charcount + 1
  end,

  Code = function(el) 
    _,n = el.text:gsub("%S+","") 
    words = words + n 
    charcount = charcount + utf8.len(el.text)
  end, 

  CodeBlock = function(el) 
    _,n = el.text:gsub("%S+","") 
    words = words + n 
    charcount = charcount + utf8.len(el.text)
  end 
} 

--[[
It looks as if the structure of the .md file is read as PARAs within BLOCKs.
All of the Blocks' character-count is contained within the Paras.

Therefore, the correct approach is to catch Blocks, then walk them to find
separately the Paras within, and count their character content separately.

At the end of each Para, I want to print its content with a character count.

At the end of each Block, I want to print the sum of its Paras' char counts.
]]

block_charcount = 0

-- Used as a 'fallback function' when no more specific function matches.
-- See https://pandoc.org/lua-filters.html#lua-filter-structure.
function Block(el)
  io.stderr:write("FOUND BLOCK: ")
  print(block_charcount .. " total characters")
end


-- Count paras in AST
function Para(el)
  io.stderr:write("  PARA: ")
  words = 0 
  charcount = 1 -- a rough approximation: assume each para adds 1 linefeed
  pandoc.walk_block(el, wordcount)
  print(charcount .. " chars")
  block_charcount = block_charcount + charcount
end

function HorizontalRule(el)
  io.stderr:write("AT HRULE: ")
  print(block_charcount .. " total characters")
  charcount_message = block_charcount .. " tweet chars"
  block_charcount = 0
  return {
    pandoc.Para(charcount_message),
    el
  }
end

