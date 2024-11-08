import { createContext, useContext, useReducer } from "react";

const PrimerDriverContext = createContext();

const initialState = {
  loadedResultsFromModule: false,
  results: {
    loaded: false,
    data: [],
  },
  mode: "",
};

const reducer = (state, action) => {
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
  const [PDState, PDDispatch] = useReducer(reducer, initialState);
  const value = { PDState, PDDispatch };

  return (
    <PrimerDriverContext.Provider value={value}>
      {props.children}
    </PrimerDriverContext.Provider>
  );
}

const PrimerDriverConsumer = PrimerDriverContext.Consumer;
const usePrimerDriverContext = () => useContext(PrimerDriverContext);
const generalContextValue = {
  PDState: initialState,
  PDDispatch: reducer,
};

export {
  PrimerDriverContext,
  PrimerDriverProvider,
  PrimerDriverConsumer,
  usePrimerDriverContext,
  generalContextValue,
};
