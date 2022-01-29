import { createContext, useReducer, useContext } from "react";

let PrimerDriverContext = createContext();

let initialState = {
  loadedResultsFromModule: false,
  results: {
    loaded: false,
    data: [],
  },
  mode: "",
};

let reducer = (state, action) => {
  switch (action.type) {
    case "updateLoadedResults": {
      return {
        ...state,
        loadedResultsFromModule: action.payload,
      };
    }
    case "updateResults": {
      return {
        ...state,
        results: action.payload,
      };
    }
    case "updateMode": {
      return {
        ...state,
        mode: action.payload,
      };
    }
    default: {
      return { ...state };
    }
  }
};

function PrimerDriverProvider(props) {
  let [PDState, PDDispatch] = useReducer(reducer, initialState);
  let value = { PDState, PDDispatch };

  return <PrimerDriverContext.Provider value={value}>{props.children}</PrimerDriverContext.Provider>;
}

let PrimerDriverConsumer = PrimerDriverContext.Consumer;
let usePrimerDriverContext = () => useContext(PrimerDriverContext);
let generalContextValue = {
  PDState: initialState,
  PDDispatch: reducer,
};

export { PrimerDriverContext, PrimerDriverProvider, PrimerDriverConsumer, usePrimerDriverContext, generalContextValue };
